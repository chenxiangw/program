
#include "compiler.h"
#include "vm.hpp"
#include <boost/foreach.hpp>
#include <boost/variant/apply_visitor.hpp>
#include <boost/assert.hpp>


namespace client { 

	namespace code_gen{

		void program::save_code(int a){
			code.push_back(a);
		}
		void program::save_value(double a){
			value.push_back(a);
		}
		void program::save_variable_name(std::string const&name) {
			variable_names.push_back(name);
		}
		void program::save_matrix(std::shared_ptr<Matrix<double>> &matrix) {
			matrix_loc.push_back(matrix);
		}

		std::pair<bool, double> const* program::find_var(std::string const& name) const{
			std::map<std::string, std::pair<bool, double>>::const_iterator i = variables.find(name);
			if (i == variables.end())return nullptr;
			else return &i->second;
		}
		void program::add_var(std::string const& name){
			variables.insert(std::pair<std::string, std::pair<bool, double>>(name, std::pair<bool, double>(false, 0)));
		}
		void program::set_var(std::string const& name, double value) {
			variables[name].first = true;
			variables[name].second = value;
		}

		std::vector<std::string> const* program::find_func(std::string const& name) const {
			std::map < std::string, std::vector<std::string> >::const_iterator i = function_calls.find(name);
			if (i == function_calls.end())return nullptr;
			else return &i->second;
		}
		std::shared_ptr<Matrix<double>> const * program::find_matrix(std::string const& name) const {
			std::map < std::string, std::shared_ptr<Matrix<double>> >::const_iterator  i = matrixes.find(name);
			if (i == matrixes.end())return nullptr;
			else return &i->second;
		}
		void program::add_matrix(std::string const& name, std::shared_ptr<Matrix<double>> &new_matrix) {
			matrixes.insert(std::pair<std::string, std::shared_ptr<Matrix<double>>>(name, new_matrix));
		}

		void program::print_assembler() const{
			std::vector<int>::const_iterator pc = code.begin();
			std::vector<double>::const_iterator value_ptr = value.begin();
			std::vector<std::string>::const_iterator variables_loc_ptr = variable_names.begin();

			while (pc != code.end()){
				switch (*pc++){
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
				case op_load_var: {
					std::cout << "op_load_var     " << std::endl;
					break;
				}
				case op_load_mat: {
					std::cout << "op_load_mat     " << std::endl;
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
				case op_store_var_plural: {
					std::cout << "op_store_var_plural   " << std::endl;
					break;
				}
				case op_store_mat: {
					std::cout << "op_store_mat   " << std::endl;
					break;
				}
				case op_dou: {
					std::cout << "op_dou      " << std::endl;
					break;
				}
				case op_mat: {
					std::cout << "op_mat      " << std::endl;
					break;
				}
				case op_param_val: {
					std::cout << "op_param_val      " << std::endl;
					break;
				}
				case op_param_mat: {
					std::cout << "op_param_mat      " << std::endl;
					break;
				}
				case op_func_call: {
					std::cout << "op_func_call      " << std::endl;
					break;
				}
				}
			}
		}

		bool compiler::operator()(double x) const{
			program.save_code(op_dou);
			program.save_value(x);
			return true;
		}

		bool compiler::operator()(ast::variable const& x) const {
			auto const* p = program.find_var(x.var);
			if (p != nullptr) {//单值变量
				if (p->first == false) {
					std::cout << "未赋值变量 " << std::endl;
					return false;
				}
				program.save_code(op_load_var);
				program.save_variable_name(x.var);
			}
			else if(program.find_matrix(x.var)){//矩阵变量
				program.save_code(op_load_mat);
				program.save_variable_name(x.var);
			}
			else {
				std::cout << "未声明变量 " << std::endl;
				return false;
			}
			return true;
		}

		bool compiler::operator()(ast::vector const& x) const {
			int n = x.values.size();
			std::shared_ptr<Matrix<double>> new_matrix;
			if (n != 0) {
				std::list<double>::const_iterator val = x.values.begin();
				new_matrix = std::shared_ptr<Matrix<double>>(new Matrix<double>(1, n));
				for (int column = 0; val != x.values.end(); ++val, ++column) {
					(*new_matrix)[0][column] = *val;
				}
			}
			else {
				new_matrix = std::shared_ptr<Matrix<double>>(new Matrix<double>());
			}
			program.save_code(op_mat);
			program.save_matrix(new_matrix);
			return true;
		}

		bool compiler::operator()(ast::matrix const& x) const {
			int m = x.rhs.size();
			std::shared_ptr<Matrix<double>> new_matrix;
			if (m != 0) {
				std::list<client::ast::vector>::const_iterator vec = x.rhs.begin();
				int n = vec->values.size();
				new_matrix = std::shared_ptr<Matrix<double>>(new Matrix<double>(m, n));
				for (int row = 0; vec != x.rhs.end(); ++vec, ++row) {
					if (vec->values.size() != n) {
						std::cout << "矩阵列数不一致" << std::endl;
						return false;
					}
					std::list<double>::const_iterator val = vec->values.begin();
					for (int column = 0; val != vec->values.end(); ++val, ++column) {
						(*new_matrix)[row][column] = *val;
					}
				}
			}
			else {
				new_matrix = std::shared_ptr<Matrix<double>>(new Matrix<double>());
			}
			program.save_code(op_mat);
			program.save_matrix(new_matrix);
			return true;
		}

		bool compiler::operator()(ast::operation const& x) const{
			if (!boost::apply_visitor(*this, x.operand_))
				return false;
			switch (x.operator_){
				case '+': program.save_code(op_add); break;
				case '-': program.save_code(op_sub); break;
				case '*': program.save_code(op_mul); break;
				case '/': program.save_code(op_div); break;
				case '^': program.save_code(op_pow); break;
				default: BOOST_ASSERT(0); return false;
			}
			return true;
		}

		bool compiler::operator()(ast::signed_ const& x) const{
			if (!boost::apply_visitor(*this, x.operand_))
				return false;
			switch (x.sign)
			{
				case '-': program.save_code(op_neg); break;
				case '+': break;
				default: BOOST_ASSERT(0); return false;
			}
			return true;
		}

		bool compiler::operator()(ast::function_call const& x)const {
			std::vector<std::string> const *param_found = program.find_func(x.function_name);
			if (param_found ==nullptr){
				std::cout << "不存在函数" << std::endl;
				return false;
			}
			std::vector<std::string>::const_iterator param_found_ptr = param_found->begin();
			if (param_found->size() != x.args.size()) {
				std::cout << "参数个数不匹配" << std::endl;
				return false;
			}
			BOOST_FOREACH(ast::parameter const& param, x.args) {
				if (typeid(param.param) == typeid(double)) {//参数为double
					if (typeid(double).name() != *param_found_ptr++) {//匹配参数不是double
						std::cout << "参数类型不匹配" << std::endl;
						return false;
					}
					program.save_code(op_param_val);
					program.save_value(boost::get<double>(param.param));
				}
				else {//参数为变量
					std::string param_type;
					auto &variable = boost::get<client::ast::variable>(param.param);
					variable.var;
					auto const* p = program.find_var(variable.var);
					if (p != nullptr) {//参数为单值变量
						if (p->first == false) {
							std::cout << "未赋值变量 " << std::endl;
							return false;
						}
						if (typeid(double).name() != *param_found_ptr++) {//匹配参数不是double
							std::cout << "参数类型不匹配" << std::endl;
							return false;
						}
						program.save_code(op_param_val);
						program.save_value(p->second);
					}
					else {
						auto const q = program.find_matrix(variable.var);
						if (q != nullptr) {//参数为矩阵变量
							if (typeid(std::shared_ptr<Matrix<double>>).name() != *param_found_ptr++) {//匹配参数不是矩阵
								std::cout << "参数类型不匹配" << std::endl;
								return false;
							}
							program.save_code(op_param_mat);
							program.save_matrix(std::shared_ptr<Matrix<double>>(*q));
						}
						else {
							std::cout << "未声明变量 " << std::endl;
							return false;
						}
					}
				}
			}
			program.save_code(op_func_call);
			program.save_variable_name(x.function_name);
			return true;
		}

		bool compiler::operator()(ast::expression const& x) const{
			if (!boost::apply_visitor(*this, x.first))
				return false;
			BOOST_FOREACH(ast::operation const& oper, x.rest){
				if (!(*this)(oper))
					return false;
			}
			return true;
		}

		bool compiler::operator()(ast::assignment const& x) const{
			if (!(*this)(x.rhs))return false;
			auto const* p = program.find_var(x.lhs);
			if (p != nullptr) {//单值变量
				if (p->first == false) {
					std::cout << "未赋值变量 " << std::endl;
					return false;
				}
				program.save_code(op_store_var);
				program.save_variable_name(x.lhs);
			}
			else if (program.find_matrix(x.lhs)) {//矩阵变量
				program.save_code(op_store_mat);
				program.save_variable_name(x.lhs);
			}
			else {//根据右值确定类型
				program.save_code(op_store);
				program.save_variable_name(x.lhs);
			}
			return true;
		}

		bool compiler::operator()(ast::variable_declaration const& x) const {
			if (x.rhs.size()>0) {//存在赋值语句
				BOOST_FOREACH(ast::expression const& expr, x.rhs){
					if (!(*this)(expr))
						return false;
				}
				int count_var = 0;
				for (auto &elem : x.vars) {
					auto const*p = program.find_var(elem);
					if (p != nullptr) {
						std::cout << "已存在变量 " << std::endl;
						continue;
					}
					++count_var;
					program.add_var(elem);
					program.save_variable_name(elem);
				}
				if (count_var > 0) {
					if (count_var > 1) {
						for (int i = 1; i < count_var; ++i) {
							program.save_code(op_store_var_plural);
						}
					}
					program.save_code(op_store_var);
				}
			}
			else {
				for (auto &elem : x.vars) {
					auto const*p = program.find_var(elem);
					if (p != nullptr) {
						std::cout << "已存在变量 " << std::endl;
						continue;
					}
					program.add_var(elem);
				}
			}
			return true;

		}

		bool compiler::operator()(ast::statement_list const& x) const{
			BOOST_FOREACH(ast::statement const& s, x){
				if (!boost::apply_visitor(*this, s)){
					program.clear();
					return false;
				}
			}
			return true;
		}
	}
}

