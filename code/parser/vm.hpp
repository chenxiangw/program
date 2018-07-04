
#pragma once

#include"compiler.h"
#include<cmath>
#include <vector>

namespace client{

    ///////////////////////////////////////////////////////////////////////////
    //  The Virtual Machine
    ///////////////////////////////////////////////////////////////////////////
    
    class vmachine{

    public:
        vmachine(Global_variable &global_variable)
        : global_variable(global_variable){
        }

		bool execute(client::program &program) {
			std::vector<int>::const_iterator code_ptr = program.get_code().begin();
			std::vector<int>::const_iterator code_end = program.get_code().end();
			std::vector<std::string>::const_iterator variable_names_ptr = program.get_variable_names().begin();
			std::vector<data_type>::const_iterator data_ptr = program.get_data().begin();

			while (code_ptr != code_end) {
				switch (*code_ptr++) {
				case op_neg: {
					if (stack.back().type() == typeid(std::shared_ptr<double>)) {
						std::shared_ptr<double> &value = boost::get<std::shared_ptr<double>>(stack.back());
						*value = -*value;
					}
					else if (stack.back().type() == typeid(std::shared_ptr<Matrix<double>>)) {
						std::shared_ptr<Matrix<double>> &matrix = boost::get<std::shared_ptr<Matrix<double>>>(stack.back());
						matrix = matrix->reverse_sign();
					}
					break;
				}
				case op_add: {
					if (stack.back().type() == typeid(std::shared_ptr<double>)) {
						std::shared_ptr<double> value2 = boost::get<std::shared_ptr<double>>(stack.back());
						stack.pop_back();
						if (stack.back().type() == typeid(std::shared_ptr<double>)) {//double+double
							std::shared_ptr<double> &value1 = boost::get<std::shared_ptr<double>>(stack.back());
							*value1 += *value2;
						}
						else if (stack.back().type() == typeid(std::shared_ptr<Matrix<double>>)) {//matrix+double
							std::shared_ptr<Matrix<double>> &matrix1 = boost::get<std::shared_ptr<Matrix<double>>>(stack.back());
							matrix1 = *matrix1 + *value2;
						}
					}
					else if (stack.back().type() == typeid(std::shared_ptr<Matrix<double>>)) {
						std::shared_ptr<Matrix<double>> matrix2 = boost::get<std::shared_ptr<Matrix<double>>>(stack.back());
						stack.pop_back();
						if (stack.back().type() == typeid(std::shared_ptr<double>)) {//double+matrix
							std::shared_ptr<double> &value1 = boost::get<std::shared_ptr<double>>(stack.back());
							matrix2 = *matrix2 + *value1;
							stack.pop_back();
							stack.push_back(matrix2);
						}
						else if (stack.back().type() == typeid(std::shared_ptr<Matrix<double>>)) {//matrix+matrix
							std::shared_ptr<Matrix<double>> &matrix1 = boost::get<std::shared_ptr<Matrix<double>>>(stack.back());
							if (matrix1->msize() != matrix2->msize() || matrix1->nsize() != matrix2->nsize()) {
								std::cout << "矩阵维度不正确" << std::endl;
								return false;
							}
							matrix1 = *matrix1 + *matrix2;
						}
					}
					break;
				}
				case op_sub: {
					if (stack.back().type() == typeid(std::shared_ptr<double>)) {
						std::shared_ptr<double> value2 = boost::get<std::shared_ptr<double>>(stack.back());
						stack.pop_back();
						if (stack.back().type() == typeid(std::shared_ptr<double>)) {//double+double
							std::shared_ptr<double> &value1 = boost::get<std::shared_ptr<double>>(stack.back());
							*value1 -= *value2;
						}
						else if (stack.back().type() == typeid(std::shared_ptr<Matrix<double>>)) {//matrix+double
							std::shared_ptr<Matrix<double>> &matrix1 = boost::get<std::shared_ptr<Matrix<double>>>(stack.back());
							matrix1 = *matrix1 - *value2;
						}
					}
					else if (stack.back().type() == typeid(std::shared_ptr<Matrix<double>>)) {
						std::shared_ptr<Matrix<double>> matrix2 = boost::get<std::shared_ptr<Matrix<double>>>(stack.back());
						stack.pop_back();
						if (stack.back().type() == typeid(std::shared_ptr<double>)) {//double+matrix
							std::shared_ptr<double> &value1 = boost::get<std::shared_ptr<double>>(stack.back());
							matrix2 = *matrix2 - *value1;
							stack.pop_back();
							stack.push_back(matrix2);
						}
						else if (stack.back().type() == typeid(std::shared_ptr<Matrix<double>>)) {//matrix+matrix
							std::shared_ptr<Matrix<double>> &matrix1 = boost::get<std::shared_ptr<Matrix<double>>>(stack.back());
							if (matrix1->msize() != matrix2->msize() || matrix1->nsize() != matrix2->nsize()) {
								std::cout << "矩阵维度不正确" << std::endl;
								return false;
							}
							matrix1 = *matrix1 - *matrix2;
						}
					}
					break;
				}
				case op_mul: {
					if (stack.back().type() == typeid(std::shared_ptr<double>)) {
						std::shared_ptr<double> value2 = boost::get<std::shared_ptr<double>>(stack.back());
						stack.pop_back();
						if (stack.back().type() == typeid(std::shared_ptr<double>)) {//double*double
							std::shared_ptr<double> &value1 = boost::get<std::shared_ptr<double>>(stack.back());
							*value1 *= *value2;
						}
						else if (stack.back().type() == typeid(std::shared_ptr<Matrix<double>>)) {//matrix*double
							std::shared_ptr<Matrix<double>> &matrix1 = boost::get<std::shared_ptr<Matrix<double>>>(stack.back());
							matrix1 = *matrix1 * *value2;
						}
					}
					else if (stack.back().type() == typeid(std::shared_ptr<Matrix<double>>)) {
						std::shared_ptr<Matrix<double>> matrix2 = boost::get<std::shared_ptr<Matrix<double>>>(stack.back());
						stack.pop_back();
						if (stack.back().type() == typeid(std::shared_ptr<double>)) {//double*matrix
							std::shared_ptr<double> &value1 = boost::get<std::shared_ptr<double>>(stack.back());
							matrix2 = *matrix2 * *value1;
							stack.pop_back();
							stack.push_back(matrix2);
						}
						else  if (stack.back().type() == typeid(std::shared_ptr<Matrix<double>>)) {//matrix*matrix
							std::shared_ptr<Matrix<double>> &matrix1 = boost::get<std::shared_ptr<Matrix<double>>>(stack.back());
							if (matrix1->nsize() != matrix2->msize()) {
								std::cout << "矩阵维度不正确" << std::endl;
								return false;
							}
							matrix1 = *matrix1 * *matrix2;
						}
					}
					break;
				}
				case op_div: {
					if (stack.back().type() == typeid(std::shared_ptr<double>)) {
						std::shared_ptr<double> value2 = boost::get<std::shared_ptr<double>>(stack.back());
						if (*value2 == 0) {
							std::cout << "除以0" << std::endl;
							return false;
						}
						stack.pop_back();
						if (stack.back().type() == typeid(std::shared_ptr<double>)) {//double/double
							std::shared_ptr<double> &value1 = boost::get<std::shared_ptr<double>>(stack.back());
							*value1 /= *value2;
						}
						else if (stack.back().type() == typeid(std::shared_ptr<Matrix<double>>)) {//matrix/double
							std::shared_ptr<Matrix<double>> &matrix1 = boost::get<std::shared_ptr<Matrix<double>>>(stack.back());
							matrix1 = *matrix1 / *value2;
						}
					}
					else  if (stack.back().type() == typeid(std::shared_ptr<Matrix<double>>)) {
						std::shared_ptr<Matrix<double>> matrix2 = boost::get<std::shared_ptr<Matrix<double>>>(stack.back());
						stack.pop_back();
						if (stack.back().type() == typeid(std::shared_ptr<double>)) {//double/matrix
							std::cout << "单值/矩阵未定义" << std::endl;
							return false;
						}
						else if (stack.back().type() == typeid(std::shared_ptr<Matrix<double>>)) {//matrix/matrix
							std::shared_ptr<Matrix<double>> &matrix1 = boost::get<std::shared_ptr<Matrix<double>>>(stack.back());
							if (matrix1->nsize() != matrix2->msize()) {
								std::cout << "矩阵维度不正确" << std::endl;
								return false;
							}
							matrix1 = *matrix1 / *matrix2;
							if (matrix1 == nullptr)return false;
						}
					}
					break;
				}
				case op_pow: {
					if (stack.back().type() == typeid(std::shared_ptr<double>)) {
						std::shared_ptr<double> value2 = boost::get<std::shared_ptr<double>>(stack.back());
						stack.pop_back();
						if (stack.back().type() == typeid(std::shared_ptr<double>)) {
							std::shared_ptr<double> &value1 = boost::get<std::shared_ptr<double>>(stack.back());
							*value1 = pow(*value1, *value2);
						}
						else {
							std::cout << "运算未定义^" << std::endl;
						}
					}
					else {
						std::cout << "运算未定义^" << std::endl;
					}
					break;
				}

				case op_load: {
					std::string variable_name = *variable_names_ptr++;
					data_type *variable = global_variable.get_variable(variable_name);
					if (variable != nullptr) {
						stack.push_back(*variable);
					}
					else {
						std::cout << "变量不存在" << std::endl;
						return false;
					}
					break;
				}
				case op_store: {
					bool success = global_variable.assign_variable(*variable_names_ptr++, stack.back());
					if (!success)return false;
					stack.pop_back();
					break;
				}
				case op_store_var: {
					double value = *boost::get<std::shared_ptr<double>>(stack.back());
					bool success = global_variable.assign_var(*variable_names_ptr++, std::shared_ptr<double>(new double(value)));
					if (!success)return false;
					stack.pop_back();
					break;
				}
				case op_data: {
					stack.push_back(*data_ptr++);
					break;
				}

				case op_param: {
					return false;
					break;
				}
				case op_func_call: {
					std::string function_name = *variable_names_ptr++;
					std::vector<data_type> params;
					//初始化参数
					while (code_ptr != code_end && *code_ptr == op_param) {
						++code_ptr;
						switch (*code_ptr++) {
						case op_load: {
							std::string variable_name = *variable_names_ptr++;
							data_type *variable = global_variable.get_variable(variable_name);
							if (variable != nullptr) {
								params.push_back(*variable);
							}
							else {
								std::cout << "变量不存在" << std::endl;
								return false;
							}
							break;
						}
						case op_data: {
							params.push_back(*data_ptr++);
							break;
						}
						}
					}
					auto  &builtin_functions = global_variable.get_builtin_functions();
					auto func = builtin_functions.find(function_name);
					if (func != builtin_functions.end()) {//调用内建函数
						if (func->second(params)) {
							stack.insert(stack.end(), params.begin(), params.end());
						}
						else return false;
					}
					else {//调用自定义函数
						std::shared_ptr<client::function_program> user_function = global_variable.get_user_function(function_name);
						if (user_function == nullptr)return false;
						if (!user_function->input(params))return false;
						if (!execute_user_function(*user_function, global_variable.get_builtin_functions()))return false;
						global_variable.store_output(*user_function);
						break;
					}
				}
				}
			}
			return true;
		}
		bool execute_user_function(client::function_program &program, std::map < std::string, bool(*)(std::vector<data_type>&)>&builtin_functions) {
			std::vector<int>::const_iterator code_ptr = program.get_code().begin();
			std::vector<int>::const_iterator code_end = program.get_code().end();
			std::vector<std::string>::const_iterator variable_names_ptr = program.get_variable_names().begin();
			std::vector<data_type>::const_iterator data_ptr = program.get_data().begin();
			std::vector<data_type> stack;
			while (code_ptr != code_end) {
				switch (*code_ptr++) {
				case op_neg: {
					if (stack.back().type() == typeid(std::shared_ptr<double>)) {
						std::shared_ptr<double> &value = boost::get<std::shared_ptr<double>>(stack.back());
						*value = -*value;
					}
					else if (stack.back().type() == typeid(std::shared_ptr<Matrix<double>>)) {
						std::shared_ptr<Matrix<double>> &matrix = boost::get<std::shared_ptr<Matrix<double>>>(stack.back());
						matrix = matrix->reverse_sign();
					}
					break;
				}
				case op_add: {
					if (stack.back().type() == typeid(std::shared_ptr<double>)) {
						std::shared_ptr<double> value2 = boost::get<std::shared_ptr<double>>(stack.back());
						stack.pop_back();
						if (stack.back().type() == typeid(std::shared_ptr<double>)) {//double+double
							std::shared_ptr<double> &value1 = boost::get<std::shared_ptr<double>>(stack.back());
							*value1 += *value2;
						}
						else if (stack.back().type() == typeid(std::shared_ptr<Matrix<double>>)) {//matrix+double
							std::shared_ptr<Matrix<double>> &matrix1 = boost::get<std::shared_ptr<Matrix<double>>>(stack.back());
							matrix1 = *matrix1 + *value2;
						}
					}
					else if (stack.back().type() == typeid(std::shared_ptr<Matrix<double>>)) {
						std::shared_ptr<Matrix<double>> matrix2 = boost::get<std::shared_ptr<Matrix<double>>>(stack.back());
						stack.pop_back();
						if (stack.back().type() == typeid(std::shared_ptr<double>)) {//double+matrix
							std::shared_ptr<double> &value1 = boost::get<std::shared_ptr<double>>(stack.back());
							matrix2 = *matrix2 + *value1;
							stack.pop_back();
							stack.push_back(matrix2);
						}
						else if (stack.back().type() == typeid(std::shared_ptr<Matrix<double>>)) {//matrix+matrix
							std::shared_ptr<Matrix<double>> &matrix1 = boost::get<std::shared_ptr<Matrix<double>>>(stack.back());
							if (matrix1->msize() != matrix2->msize() || matrix1->nsize() != matrix2->nsize()) {
								std::cout << "矩阵维度不正确" << std::endl;
								return false;
							}
							matrix1 = *matrix1 + *matrix2;
						}
					}
					break;
				}
				case op_sub: {
					if (stack.back().type() == typeid(std::shared_ptr<double>)) {
						std::shared_ptr<double> value2 = boost::get<std::shared_ptr<double>>(stack.back());
						stack.pop_back();
						if (stack.back().type() == typeid(std::shared_ptr<double>)) {//double+double
							std::shared_ptr<double> &value1 = boost::get<std::shared_ptr<double>>(stack.back());
							*value1 -= *value2;
						}
						else if (stack.back().type() == typeid(std::shared_ptr<Matrix<double>>)) {//matrix+double
							std::shared_ptr<Matrix<double>> &matrix1 = boost::get<std::shared_ptr<Matrix<double>>>(stack.back());
							matrix1 = *matrix1 - *value2;
						}
					}
					else if (stack.back().type() == typeid(std::shared_ptr<Matrix<double>>)) {
						std::shared_ptr<Matrix<double>> matrix2 = boost::get<std::shared_ptr<Matrix<double>>>(stack.back());
						stack.pop_back();
						if (stack.back().type() == typeid(std::shared_ptr<double>)) {//double+matrix
							std::shared_ptr<double> &value1 = boost::get<std::shared_ptr<double>>(stack.back());
							matrix2 = *matrix2 - *value1;
							stack.pop_back();
							stack.push_back(matrix2);
						}
						else if (stack.back().type() == typeid(std::shared_ptr<Matrix<double>>)) {//matrix+matrix
							std::shared_ptr<Matrix<double>> &matrix1 = boost::get<std::shared_ptr<Matrix<double>>>(stack.back());
							if (matrix1->msize() != matrix2->msize() || matrix1->nsize() != matrix2->nsize()) {
								std::cout << "矩阵维度不正确" << std::endl;
								return false;
							}
							matrix1 = *matrix1 - *matrix2;
						}
					}
					break;
				}
				case op_mul: {
					if (stack.back().type() == typeid(std::shared_ptr<double>)) {
						std::shared_ptr<double> value2 = boost::get<std::shared_ptr<double>>(stack.back());
						stack.pop_back();
						if (stack.back().type() == typeid(std::shared_ptr<double>)) {//double*double
							std::shared_ptr<double> &value1 = boost::get<std::shared_ptr<double>>(stack.back());
							*value1 *= *value2;
						}
						else if (stack.back().type() == typeid(std::shared_ptr<Matrix<double>>)) {//matrix*double
							std::shared_ptr<Matrix<double>> &matrix1 = boost::get<std::shared_ptr<Matrix<double>>>(stack.back());
							matrix1 = *matrix1 * *value2;
						}
					}
					else if (stack.back().type() == typeid(std::shared_ptr<Matrix<double>>)) {
						std::shared_ptr<Matrix<double>> matrix2 = boost::get<std::shared_ptr<Matrix<double>>>(stack.back());
						stack.pop_back();
						if (stack.back().type() == typeid(std::shared_ptr<double>)) {//double*matrix
							std::shared_ptr<double> &value1 = boost::get<std::shared_ptr<double>>(stack.back());
							matrix2 = *matrix2 * *value1;
							stack.pop_back();
							stack.push_back(matrix2);
						}
						else  if (stack.back().type() == typeid(std::shared_ptr<Matrix<double>>)) {//matrix*matrix
							std::shared_ptr<Matrix<double>> &matrix1 = boost::get<std::shared_ptr<Matrix<double>>>(stack.back());
							if (matrix1->nsize() != matrix2->msize()) {
								std::cout << "矩阵维度不正确" << std::endl;
								return false;
							}
							matrix1 = *matrix1 * *matrix2;
						}
					}
					break;
				}
				case op_div: {
					if (stack.back().type() == typeid(std::shared_ptr<double>)) {
						std::shared_ptr<double> value2 = boost::get<std::shared_ptr<double>>(stack.back());
						if (*value2 == 0) {
							std::cout << "除以0" << std::endl;
							return false;
						}
						stack.pop_back();
						if (stack.back().type() == typeid(std::shared_ptr<double>)) {//double/double
							std::shared_ptr<double> &value1 = boost::get<std::shared_ptr<double>>(stack.back());
							*value1 /= *value2;
						}
						else if (stack.back().type() == typeid(std::shared_ptr<Matrix<double>>)) {//matrix/double
							std::shared_ptr<Matrix<double>> &matrix1 = boost::get<std::shared_ptr<Matrix<double>>>(stack.back());
							matrix1 = *matrix1 / *value2;
						}
					}
					else  if (stack.back().type() == typeid(std::shared_ptr<Matrix<double>>)) {
						std::shared_ptr<Matrix<double>> matrix2 = boost::get<std::shared_ptr<Matrix<double>>>(stack.back());
						stack.pop_back();
						if (stack.back().type() == typeid(std::shared_ptr<double>)) {//double/matrix
							std::cout << "单值/矩阵未定义" << std::endl;
							return false;
						}
						else if (stack.back().type() == typeid(std::shared_ptr<Matrix<double>>)) {//matrix/matrix
							std::shared_ptr<Matrix<double>> &matrix1 = boost::get<std::shared_ptr<Matrix<double>>>(stack.back());
							if (matrix1->nsize() != matrix2->msize()) {
								std::cout << "矩阵维度不正确" << std::endl;
								return false;
							}
							matrix1 = *matrix1 / *matrix2;
							if (matrix1 == nullptr)return false;
						}
					}
					break;
				}
				case op_pow: {
					if (stack.back().type() == typeid(std::shared_ptr<double>)) {
						std::shared_ptr<double> value2 = boost::get<std::shared_ptr<double>>(stack.back());
						stack.pop_back();
						if (stack.back().type() == typeid(std::shared_ptr<double>)) {
							std::shared_ptr<double> &value1 = boost::get<std::shared_ptr<double>>(stack.back());
							*value1 = pow(*value1, *value2);
						}
						else {
							std::cout << "运算未定义^" << std::endl;
						}
					}
					else {
						std::cout << "运算未定义^" << std::endl;
					}
					break;
				}

				case op_load: {
					std::string variable_name = *variable_names_ptr++;
					data_type *variable = program.get_local_variable(variable_name);
					if (variable != nullptr) {
						stack.push_back(*variable);
					}
					else {
						std::cout << "变量不存在" << std::endl;
						return false;
					}
					break;
				}
				case op_store: {
					bool success = program.assign_local_variable(*variable_names_ptr++, stack.back());
					if (!success)return false;
					stack.pop_back();
					break;
				}
				case op_store_var: {
					double value = *boost::get<std::shared_ptr<double>>(stack.back());
					bool success = program.assign_var(*variable_names_ptr++, std::shared_ptr<double>(new double(value)));
					if (!success)return false;
					stack.pop_back();
					break;
				}
				case op_data: {
					stack.push_back(*data_ptr++);
					break;
				}

				case op_param: {
					return false;
					break;
				}
				case op_func_call: {
					std::string function_name = *variable_names_ptr++;
					std::vector<data_type> params;
					//初始化参数
					while (code_ptr != code_end && *code_ptr == op_param) {
						++code_ptr;
						switch (*code_ptr++) {
						case op_load: {
							std::string variable_name = *variable_names_ptr++;
							data_type *variable = program.get_local_variable(variable_name);
							if (variable != nullptr) {
								params.push_back(*variable);
							}
							else {
								std::cout << "变量不存在" << std::endl;
								return false;
							}
							break;
						}
						case op_data: {
							params.push_back(*data_ptr++);
							break;
						}
						}
					}
					auto func = builtin_functions.find(function_name);
					if (func != builtin_functions.end()) {//调用内建函数
						if (func->second(params)) {
							stack.insert(stack.end(), params.begin(), params.end());
						}
						else return false;
					}
					else {
						std::shared_ptr<client::function_program> user_function = global_variable.get_user_function(function_name);
						if (user_function == nullptr)return false;
						if (!user_function->input(params))return false;
						if (!execute_user_function(*user_function, builtin_functions))return false;
						global_variable.store_output(*user_function);
						break;
					}
				}
				}
			}
			return true;
		}

		void clear() { stack.clear(); }
		void print_stack() {
			std::vector<data_type>::const_iterator result = stack.begin();
			while (result != stack.end()) {
				if ((*result).type() == typeid(std::shared_ptr<double>)) {
					std::shared_ptr<double> &value = boost::get<std::shared_ptr<double>>(stack.back());
					std::cout << *value << std::endl;
				}
				else {
					std::shared_ptr<Matrix<double>> &matrix = boost::get<std::shared_ptr<Matrix<double>>>(stack.back());
					if (matrix != nullptr) {
						matrix->print();
					}
				}
				++result;
			}
		}

    private:
		std::vector<data_type> stack;
		Global_variable &global_variable;
    };

}


