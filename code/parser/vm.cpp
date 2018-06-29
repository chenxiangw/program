
#include "vm.hpp"
#include<cmath>



namespace client{

    bool vmachine::execute(client::code_gen::program &program) {
		std::vector<int>::const_iterator code_ptr = program.get_code().begin();
		std::vector<int>::const_iterator code_end = program.get_code().end();
		std::vector<double>::const_iterator value_ptr = program.get_value().begin();
		std::vector<std::string>::const_iterator variable_names_ptr = program.get_variables_loc().begin();
		std::vector< std::shared_ptr<Matrix<double>>>::const_iterator matrix_loc_ptr = program.get_matrix_loc().begin();

		while (code_ptr != code_end) {
            switch (*code_ptr++){
			case op_neg: {
				if (data.back().type() == typeid(double)) {
					double &value = boost::get<double>(data.back());
					value = -value;
				}
				else {
					std::shared_ptr<Matrix<double>> &matrix = boost::get<std::shared_ptr<Matrix<double>>>(data.back());
					matrix = matrix->reverse_sign();
				}
				break;
			}
			case op_add: {
				if (data.back().type() == typeid(double)) {
					double value2 = boost::get<double>(data.back());
					data.pop_back();
					if (data.back().type() == typeid(double)) {//double+double
						double &value1 = boost::get<double>(data.back());
						value1 += value2;
					}
					else {//matrix+double
						std::shared_ptr<Matrix<double>> &matrix1 = boost::get<std::shared_ptr<Matrix<double>>>(data.back());
						matrix1 = *matrix1 + value2;
					}
				}
				else {
					std::shared_ptr<Matrix<double>> matrix2 = boost::get<std::shared_ptr<Matrix<double>>>(data.back());
					data.pop_back();
					if (data.back().type() == typeid(double)) {//double+matrix
						double value1= boost::get<double>(data.back());
						matrix2 = *matrix2 + value1;
						data.pop_back();
						data.push_back(matrix2);
					}
					else {//matrix+matrix
						std::shared_ptr<Matrix<double>> &matrix1 = boost::get<std::shared_ptr<Matrix<double>>>(data.back());
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
				if (data.back().type() == typeid(double)) {
					double value2 = boost::get<double>(data.back());
					data.pop_back();
					if (data.back().type() == typeid(double)) {//double-double
						double &value1 = boost::get<double>(data.back());
						value1 -= value2;
					}
					else {//matrix-double
						std::shared_ptr<Matrix<double>> &matrix1 = boost::get<std::shared_ptr<Matrix<double>>>(data.back());
						matrix1 = *matrix1 - value2;
					}
				}
				else {
					std::shared_ptr<Matrix<double>> matrix2 = boost::get<std::shared_ptr<Matrix<double>>>(data.back());
					data.pop_back();
					if (data.back().type() == typeid(double)) {//double-matrix
						double value1 = boost::get<double>(data.back());
						matrix2 = *matrix2 - value1;
						data.pop_back();
						data.push_back(matrix2);
					}
					else {//matrix-matrix
						std::shared_ptr<Matrix<double>> &matrix1 = boost::get<std::shared_ptr<Matrix<double>>>(data.back());
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
				if (data.back().type() == typeid(double)) {
					double value2 = boost::get<double>(data.back());
					data.pop_back();
					if (data.back().type() == typeid(double)) {//double*double
						double &value1 = boost::get<double>(data.back());
						value1 *= value2;
					}
					else {//matrix*double
						std::shared_ptr<Matrix<double>> &matrix1 = boost::get<std::shared_ptr<Matrix<double>>>(data.back());
						matrix1 = *matrix1 * value2;
					}
				}
				else {
					std::shared_ptr<Matrix<double>> matrix2 = boost::get<std::shared_ptr<Matrix<double>>>(data.back());
					data.pop_back();
					if (data.back().type() == typeid(double)) {//double*matrix
						double value1 = boost::get<double>(data.back());
						matrix2 = *matrix2 * value1;
						data.pop_back();
						data.push_back(matrix2);
					}
					else {//matrix*matrix
						std::shared_ptr<Matrix<double>> &matrix1 = boost::get<std::shared_ptr<Matrix<double>>>(data.back());
						if (matrix1->nsize() != matrix2->msize() ) {
							std::cout << "矩阵维度不正确" << std::endl;
							return false;
						}
						matrix1 = *matrix1 * *matrix2;
					}
				}
				break;
			}
			case op_div: {
				if (data.back().type() == typeid(double)) {
					double value2 = boost::get<double>(data.back());
					if (value2 == 0) {
						std::cout << "除以0" << std::endl;
						return false;
					}
					data.pop_back();
					if (data.back().type() == typeid(double)) {//double/double
						double &value1 = boost::get<double>(data.back());
						value1 /= value2;
					}
					else {//matrix/double
						std::shared_ptr<Matrix<double>> &matrix1 = boost::get<std::shared_ptr<Matrix<double>>>(data.back());
						matrix1 = *matrix1 / value2;
					}
				}
				else {
					std::shared_ptr<Matrix<double>> matrix2 = boost::get<std::shared_ptr<Matrix<double>>>(data.back());
					data.pop_back();
					if (data.back().type() == typeid(double)) {//double/matrix
						std::cout << "单值/矩阵未定义" << std::endl;
						return false;
					}
					else {//matrix/matrix
						std::shared_ptr<Matrix<double>> &matrix1 = boost::get<std::shared_ptr<Matrix<double>>>(data.back());
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
				if (data.back().type() == typeid(double)) {
					double value2 = boost::get<double>(data.back());
					data.pop_back();
					if (data.back().type() == typeid(double)) {
						double &value1 = boost::get<double>(data.back());
						value1 = pow(value1, value2);
					}
					else {
						std::cout << "矩阵运算未定义^" << std::endl;
					}
				}
				else {
					std::cout << "矩阵运算未定义^" << std::endl;
				}
				break;
			}
			
			case op_load_var: {
				data.push_back(variables[*variable_names_ptr++].second);
				break;
			}
			case op_load_mat: {
				data.push_back(matrixes[*variable_names_ptr++]);
				break;
			}
			case op_store: {
				if (data.back().type() == typeid(double)) {
					variables.insert(
						std::pair<std::string, std::pair<bool, double>>(*variable_names_ptr++, 
							std::pair<bool, double>(true, boost::get<double>(data.back()))));
					data.pop_back();
				}
				else if(data.back().type() == typeid(std::shared_ptr<Matrix<double>>)){
					matrixes.insert(std::pair< std::string, std::shared_ptr<Matrix<double>>>(*variable_names_ptr++, boost::get<std::shared_ptr<Matrix<double>>>(data.back())));
					data.pop_back();
				}
				break;
			}
			case op_store_var: {
				if (data.back().type() != typeid(double))return false;
				variables[*variable_names_ptr++] = std::pair<bool, double>(true, boost::get<double>(data.back()));
				data.pop_back();
				break;
			}
			case op_store_var_plural: {
				if (data.back().type() != typeid(double))return false;
				variables[*variable_names_ptr++] = std::pair<bool, double>(true, boost::get<double>(data.back()));
				break;
			}
			case op_store_mat: {
				if (data.back().type() != typeid(std::shared_ptr<Matrix<double>>))return false;
				std::string matrix_name = *variable_names_ptr++;
				if (matrixes.find(matrix_name) == matrixes.end()) {
					matrixes.insert(std::pair< std::string, std::shared_ptr<Matrix<double>>>(matrix_name, boost::get<std::shared_ptr<Matrix<double>>>(data.back())));
				}
				else {
					matrixes[matrix_name] = boost::get<std::shared_ptr<Matrix<double>>>(data.back());
				}
				data.pop_back();
				break;
			}

			case op_dou: {
				data.push_back(*value_ptr++);
				break;
			}
			case op_mat: {
				data.push_back(*matrix_loc_ptr++);
				break;
			}

			case op_param_val: {
				data.push_back(*value_ptr++);
				break;
			}
			case op_param_mat: {
				data.push_back(*matrix_loc_ptr++);
				break;
			}
			case op_func_call: {
				builtin_functions[*variable_names_ptr++](this);
				break;
			}
            }
        }
		return true;
    }


	void vmachine::print_stack() {
		std::vector<boost::variant<double,
			std::shared_ptr<Matrix<double>>
			>>::const_iterator result = data.begin();
		while (result != data.end()) {
			if ((*result).type() == typeid(double)) {
				double &value = boost::get<double>(data.back());
				std::cout << value << std::endl;
			}
			else {
				std::shared_ptr<Matrix<double>> &matrix = boost::get<std::shared_ptr<Matrix<double>>>(data.back());
				if (matrix != nullptr) {
					matrix->print();
				}
			}
			++result;
		}
	}
}


