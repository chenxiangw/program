#pragma once
#include<set>
#include"builtin_functions.h"


class Add_builtin_function {
public:
	Add_builtin_function(std::set<std::string> &illegal_names, 
		std::map < std::string, bool(*)(std::vector<data_type>&data)> &builtin_functions)
		:illegal_names(illegal_names), builtin_functions(builtin_functions){
	}
	void add(std::string name, bool(*address)(std::vector<data_type>&data)) {
		builtin_functions.insert(std::pair < std::string, bool(*)(std::vector<data_type>&data)>(name, address));
		illegal_names.insert(name);
	}
private:
	std::map < std::string, bool(*)(std::vector<data_type>&data)> &builtin_functions;
	std::set<std::string> &illegal_names;
};

class Global_variable {
public:
	void initiate() {
		illegal_names.insert({ "var","matrix","function" });
		class Add_builtin_function add_builtin_function(illegal_names, builtin_functions);
		std::string double_type = "double";
		std::string matrix_type = "matrix";
		std::string function_type = "function";

		add_builtin_function.add("cos",  cos_builtin);
		add_builtin_function.add("sin",  sin_builtin);
		add_builtin_function.add("tan",  tan_builtin);
		add_builtin_function.add("acos",  acos_builtin);
		add_builtin_function.add("asin",  asin_builtin);
		add_builtin_function.add("atan",  atan_builtin);

		add_builtin_function.add("gamma",  gamma_builtin);
		add_builtin_function.add("ln_gamma",  ln_gamma_builtin);
		add_builtin_function.add("gamma_P",  gamma_P_builtin);
		add_builtin_function.add("gamma_Q",  gamma_Q_builtin);
		add_builtin_function.add("factorial",  factorial_builtin);
		add_builtin_function.add("binomial", binomial_builtin);
		add_builtin_function.add("beta", beta_builtin);
		add_builtin_function.add("incomplete_beta", incomplete_beta_builtin);
		add_builtin_function.add("error_function",  error_function_builtin);
		add_builtin_function.add("Legendre_p", Legendre_p_builtin);
		add_builtin_function.add("Legendre_q", Legendre_q_builtin);

		add_builtin_function.add("inverse", inverse_builtin);
		add_builtin_function.add("determinant", determinant_builtin);
		add_builtin_function.add("solve", solve_builtin);
		illegal_names_sys = illegal_names;
	}

	bool declara_var(std::string const &name) {
		if (illegal_names.find(name) != illegal_names.end()) {
			std::cout << "变量名已存在或为关键字" << std::endl;
			return false;
		}
		illegal_names.insert(name);
		variables.insert(std::pair<std::string, data_type>(name, std::shared_ptr<double>(nullptr)));
		return true;
	}

	bool assign_var(std::string const &name, std::shared_ptr<double> &var) {
		auto variable = (data_type)var;
		if (illegal_names.find(name) == illegal_names.end()) {
			illegal_names.insert(name);
			variables.insert(std::pair<std::string, data_type>(name, variable));
		}
		else {
			auto result = variables.find(name);
			if (result == variables.end()) {
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
	bool assign_function(std::string const &name, std::shared_ptr<client::function_program> &function) {
		auto variable = (data_type)function;
		if (illegal_names.find(name) == illegal_names.end()) {
			illegal_names.insert(name);
			variables.insert(std::pair<std::string, data_type>(name, variable));
		}
		else {
			auto result = variables.find(name);
			if (result == variables.end()) {
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
	bool assign_variable(std::string const &name, data_type &variable) {
		if (illegal_names.find(name) == illegal_names.end()) {
			illegal_names.insert(name);
			variables.insert(std::pair<std::string, data_type>(name, variable));
		}
		else{
			auto result = variables.find(name);
			if (result == variables.end()) {
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

	std::set < std::string > const &get_illegal_names_sys()const { return illegal_names_sys; }
	data_type *get_variable(std::string const &name) {
		auto result = variables.find(name);
		if (result == variables.end())return nullptr;
		return &(result->second);
	}
	std::map < std::string, bool(*)(std::vector<data_type>&)>& get_builtin_functions() {
		return builtin_functions;
	}
	std::shared_ptr<client::function_program> get_user_function(std::string const &name) {
		if (variables.find(name) == variables.end() || variables[name].type() != typeid(std::shared_ptr<client::function_program>)) {
			std::cout << "不存在函数" << std::endl;
			return nullptr;
		}
		else {
			return boost::get<std::shared_ptr<client::function_program>>(variables[name]);
		}
	}

	bool store_output(client::function_program &user_function) {
		std::vector<std::pair<std::string, data_type>> output;
		user_function.get_output(output);
		for (auto &elem : output) {
			if (!assign_variable(elem.first, elem.second)) {
				return false;
			}
		}
		return true;
	}

private:
	std::map < std::string, data_type> variables;
	std::map < std::string, bool(*)(std::vector<data_type>&)> builtin_functions;
	std::set < std::string > illegal_names;
	std::set < std::string > illegal_names_sys;
};
