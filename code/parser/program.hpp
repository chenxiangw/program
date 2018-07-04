#pragma once
#include "../core/Matrix.hpp"
#include <boost/variant/recursive_variant.hpp>
#include <vector>
#include<set>
#include <map>
#include <string>
namespace client {
	struct function_program;
}
typedef boost::variant<
	std::shared_ptr<double>,
	std::shared_ptr<Matrix<double>>,
	std::shared_ptr<client::function_program>
> data_type;

namespace client {

	enum byte_code {
		op_neg,
		op_add,
		op_sub,
		op_mul,
		op_div,
		op_pow,

		op_load,    //  load a variable
		op_store,
		op_store_var,   //  store a var
		op_data,

		op_param,
		op_func_call,
	};

	///////////////////////////////////////////////////////////////////////////
	//  The Program
	///////////////////////////////////////////////////////////////////////////
	struct program {

		void save_code(int a) {
			code.push_back(a);
		}
		void save_data(double x) {
			data.push_back(std::shared_ptr<double>(new double(x)));
		}
		void save_data(std::shared_ptr<Matrix<double>> &matrix) {
			data.push_back(matrix);
		}
		void save_variable_name(std::string const&name) {
			variable_names.push_back(name);
		}

		std::vector<int> const& get_code() const { return code; }
		std::vector<data_type> const &get_data() const{ return data; }
		std::vector<std::string> const& get_variable_names() const { return variable_names; }

		void clear() {
			code.clear();
			data.clear();
			variable_names.clear();
		}
		void print_assembler() const {
			std::vector<int>::const_iterator code_ptr = code.begin();
			std::vector<std::string>::const_iterator variable_names_ptr = variable_names.begin();
			std::vector<data_type>::const_iterator data_ptr = data.begin();
			while (code_ptr != code.end()) {
				switch (*code_ptr++) {
				case op_neg: {
					std::cout << "op_neg" << std::endl;
					break;
				}
				case op_add: {
					std::cout << "op_add" << std::endl;
					break;
				}
				case op_sub: {
					std::cout << "op_sub" << std::endl;
					break;
				}
				case op_mul: {
					std::cout << "op_mul" << std::endl;
					break;
				}
				case op_div: {
					std::cout << "op_div" << std::endl;
					break;
				}
				case op_pow: {
					std::cout << "op_pow" << std::endl;
					break;
				}
				case op_load: {
					std::cout << "op_load     " << std::endl;
					break;
				}
				case op_store: {
					std::cout << "op_store    " << std::endl;
					break;
				}
				case op_store_var: {
					std::cout << "op_store_var    " << std::endl;
					break;
				}
				case op_data: {
					std::cout << "op_data      " << std::endl;
					break;
				}
				case op_param: {
					std::cout << "op_param      " << std::endl;
					break;
				}
				case op_func_call: {
					std::cout << "op_func_call      " << std::endl;
					break;
				}
				}
			}
		}

	private:

		std::vector<int> code;
		std::vector<data_type> data;
		std::vector<std::string> variable_names;
	};

	///////////////////////////////////////////////////////////////////////////
	//  The function_program
	///////////////////////////////////////////////////////////////////////////
	struct function_program {


		void save_code(int a) {
			code.push_back(a);
		}
		void save_data(double x) {
			data.push_back(std::shared_ptr<double>(new double(x)));
		}
		void save_data(std::shared_ptr<Matrix<double>> &matrix) {
			data.push_back(matrix);
		}
		void save_variable_name(std::string const&name) {
			variable_names.push_back(name);
		}

		bool declara_var(std::string const &name) {
			if (illegal_names.find(name) != illegal_names.end()) {
				std::cout << "变量名已存在或为关键字" << std::endl;
				return false;
			}
			illegal_names.insert(name);
			local_variables.insert(std::pair<std::string, data_type>(name, std::shared_ptr<double>(nullptr)));
			return true;
		}
		bool add_var(std::string const& name) {
			if (illegal_names.find(name) != illegal_names.end())return false;
			illegal_names.insert(name);
			//local_variable.insert(std::pair<std::string, data_type >(name, std::pair<bool, double>(false, 0)));
			return true;
		}
		bool add_matrix(std::string const& name) {
			if (illegal_names.find(name) != illegal_names.end())return false;
			illegal_names.insert(name);
			//local_variable.insert(std::pair<std::string, data_type >(name, std::shared_ptr<Matrix<double>>(nullptr)));
			return true;
		}

		bool assign_var(std::string const &name, std::shared_ptr<double> &var) {
			auto variable = (data_type)var;
			if (illegal_names.find(name) == illegal_names.end()) {
				illegal_names.insert(name);
				local_variables.insert(std::pair<std::string, data_type>(name, variable));
			}
			else {
				auto result = local_variables.find(name);
				if (result == local_variables.end()) {
					std::cout << "变量名为关键字" << std::endl;
					return false;
				}
				if (result->second.type() != variable.type()) {
					std::cout << "赋值类型不正确" << std::endl;
					return false;
				}
				result->second = variable;
			}
			return true;
		}
		bool assign_local_variable(std::string const &name, data_type &variable) {
			if (illegal_names.find(name) == illegal_names.end()) {
				illegal_names.insert(name);
				local_variables.insert(std::pair<std::string, data_type>(name, variable));
			}
			else {
				auto result = local_variables.find(name);
				if (result == local_variables.end()) {
					std::cout << "变量名为关键字" << std::endl;
					return false;
				}
				if (result->second.type() != variable.type()) {
					std::cout << "赋值类型不正确" << std::endl;
					return false;
				}
				result->second = variable;
			}
			return true;
		}

		data_type *get_local_variable(std::string const &name) {
			auto result = local_variables.find(name);
			if (result == local_variables.end())return nullptr;
			return &(result->second);
		}
		bool get_output(std::vector<std::pair<std::string,data_type>> &output) {
			for (auto &elem : output_name) {
				if (local_variables.find(elem) == local_variables.end()) {
					std::cout << "函数输出变量未定义" << std::endl;
					return false;
				}
				output.push_back(std::pair<std::string, data_type>(elem, local_variables[elem]));
			}
			return true;
		}
		std::vector<int> const& get_code() const { return code; }
		std::vector<data_type> const &get_data() const { return data; }
		std::vector<std::string> const& get_variable_names() const { return variable_names; }

		bool initiate(std::string const &content_, std::vector<std::string> const &output_name_, std::vector<std::pair<std::string, std::string>> const&input_type_, std::set < std::string >const &illegal_names_sys) {
			content = content_;
			illegal_names = illegal_names_sys;
			for (auto &elem : input_type_) {
				if (illegal_names.find(elem.first) == illegal_names.end()) {
					illegal_names.insert(elem.first);
					if (elem.second == typeid(std::shared_ptr<double>).name()) {
						add_var(elem.first);
					}
					else if (elem.second == typeid(std::shared_ptr<Matrix<double>>).name()) {
						add_matrix(elem.first);
					}
					else { return false; }
				}
				else {
					return false;
				}
			}
			for (auto &elem : output_name_) {
				if (illegal_names.find(elem) != illegal_names.end())return false;
			}
			output_name = output_name_;
			input_type = input_type_;
			return true;
		}

		bool input(std::vector<data_type> &params) {
			if (params.size() != input_type.size()) {
				std::cout << "函数参数个数不匹配" << std::endl;
				return false;
			}
			for (size_t i = 0; i < params.size(); ++i) {
				if (params[i].type().name() != input_type[i].second) {
					std::cout << "函数参数类型不匹配" << std::endl;
					return false;
				}
				local_variables[input_type[i].first] = params[i];
			}
			return true;
		}

	private:
		std::string content;
		std::vector<std::string> output_name;
		std::vector<std::pair<std::string, std::string>> input_type;

		std::map<std::string, data_type> local_variables;
		std::set < std::string > illegal_names;

		std::vector<int> code;
		std::vector<data_type> data;
		std::vector<std::string> variable_names;
	};
}