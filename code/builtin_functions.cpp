#include<iostream>
#include<string>
#include<map>
#include<utility>
#include<vector>
#include<limits>
#include<cmath>

#include "builtin_functions.h"
#include "core/special_functions.h"
//#include"Evaluation_function.h"
//#include"./core/numerical_recipes.h"


void inverse_builtin(client::vmachine  *vm) {
	std::vector<
		boost::variant<double, std::shared_ptr<Matrix<double>>>
	>&data = vm->get_data();
	if (data.back().type() != typeid(std::shared_ptr<Matrix<double>>))return;
	std::shared_ptr<Matrix<double>> &matrix = boost::get< std::shared_ptr<Matrix<double>>>(data.back());
	matrix = matrix->inverse();
	return;
}
void determinant_builtin(client::vmachine  *vm){
	std::vector<
		boost::variant<double, std::shared_ptr<Matrix<double>>>
	>&data = vm->get_data();
	if (data.back().type() != typeid(std::shared_ptr<Matrix<double>>))return;
	std::shared_ptr<Matrix<double>> &matrix = boost::get< std::shared_ptr<Matrix<double>>>(data.back());
	double result = matrix->determinant();
	data.pop_back();
	data.push_back(result);
}
void solve_builtin(client::vmachine  *vm) {
	std::vector<
		boost::variant<double, std::shared_ptr<Matrix<double>>>
	>&data = vm->get_data();
	if (data.back().type() != typeid(std::shared_ptr<Matrix<double>>))return;
	std::shared_ptr<Matrix<double>> b = boost::get< std::shared_ptr<Matrix<double>>>(data.back());
	data.pop_back();
	std::shared_ptr<Matrix<double>> A = boost::get< std::shared_ptr<Matrix<double>>>(data.back());
	data.pop_back();
	if (A->msize() != b->msize()) {
		std::cout << "矩阵格式错误" << std::endl;
		return;
	}
	auto result = solve_linear_SVD_method(*A, *b);
	int question_num = result.size();
	for (int i = 0; i < question_num; ++i) {
		std::cout << "第" << i + 1 << "个方程：" << std::endl;
		result[i]->print();
	}
}

//void integrate(const std::vector<std::string> *input) {
//	using namespace std;
//	if (input->size() != 3) {
//		cout << "参数错误" << endl;
//		return;
//	}
//	string function_name = (*input)[0];
//	if (functions.find(function_name) == functions.end()) {
//		string error = "不存在函数" + function_name;
//		throw  error;
//	}
//	double lower_bound, upper_bound;
//	if ((*input)[1] == "-inf") lower_bound = -DBL_MAX;
//	else lower_bound = stod((*input)[1]);
//	if ((*input)[2] == "inf") upper_bound = DBL_MAX;
//	else upper_bound = stod((*input)[2]);
//	Evaluation_function func(function_name);
//
//	double result = mid_trapezoidal_integration_iter_Romberg(func, lower_bound, upper_bound);
//	function_params.push_back(to_string(result));
//	return;
//}
//void diff(const std::vector<std::string> *input) {
//	using namespace std;
//	if (input->size() != 2) {
//		cout << "参数错误" << endl;
//		return;
//	}
//	string function_name = (*input)[0];
//	if (functions.find(function_name) == functions.end()) {
//		string error = "不存在函数" + function_name;
//		throw  error;
//	}
//	double x = stod((*input)[1]);
//	Evaluation_function func(function_name);
//	double result;
//	double error = derivation_Ridders(func, result, x, 0.1);
//	function_params.push_back(to_string(result));
//	return;
//}

void cos_builtin(client::vmachine  *vm) {
	std::vector<
		boost::variant<double, std::shared_ptr<Matrix<double>>>
	>&data = vm->get_data();
	if (data.back().type() != typeid(double))return;
	double &x = boost::get<double>(data.back());
	x = cos(x);
	return;
}
void sin_builtin(client::vmachine  *vm) {
	std::vector<
		boost::variant<double, std::shared_ptr<Matrix<double>>>
	>&data = vm->get_data();
	if (data.back().type() != typeid(double))return;
	double &x = boost::get<double>(data.back());
	x = sin(x);
	return;
}
void tan_builtin(client::vmachine  *vm) {
	std::vector<
		boost::variant<double, std::shared_ptr<Matrix<double>>>
	>&data = vm->get_data();
	if (data.back().type() != typeid(double))return;
	double &x = boost::get<double>(data.back());
	x = tan(x);
	return;
}
void acos_builtin(client::vmachine  *vm) {
	std::vector<
		boost::variant<double, std::shared_ptr<Matrix<double>>>
	>&data = vm->get_data();
	if (data.back().type() != typeid(double))return;
	double &x = boost::get<double>(data.back());
	x = acos(x);
	return;
}
void asin_builtin(client::vmachine  *vm) {
	std::vector<
		boost::variant<double, std::shared_ptr<Matrix<double>>>
	>&data = vm->get_data();
	if (data.back().type() != typeid(double))return;
	double &x = boost::get<double>(data.back());
	x = asin(x);
	return;
}
void atan_builtin(client::vmachine  *vm) {
	std::vector<
		boost::variant<double, std::shared_ptr<Matrix<double>>>
	>&data = vm->get_data();
	if (data.back().type() != typeid(double))return;
	double &x = boost::get<double>(data.back());
	x = atan(x);
	return;
}

void gamma_builtin(client::vmachine  *vm) {
	std::vector<
		boost::variant<double, std::shared_ptr<Matrix<double>>>
	>&data = vm->get_data();
	if (data.back().type() != typeid(double))return;
	double &x = boost::get<double>(data.back());
	x = gamma(x);
	return;
}
void ln_gamma_builtin(client::vmachine  *vm) {
	std::vector<
		boost::variant<double, std::shared_ptr<Matrix<double>>>
	>&data = vm->get_data();
	if (data.back().type() != typeid(double))return;
	double &x = boost::get<double>(data.back());
	x = ln_gamma(x);
	return;
}
void gamma_P_builtin(client::vmachine  *vm) {
	std::vector<
		boost::variant<double, std::shared_ptr<Matrix<double>>>
	>&data = vm->get_data();
	if (data.back().type() != typeid(double))return;
	double y = boost::get<double>(data.back());
	data.pop_back();
	if (data.back().type() != typeid(double))return;
	double &x = boost::get<double>(data.back());
	x = gamma_P(x, y);
	return;
}
void gamma_Q_builtin(client::vmachine  *vm) {
	std::vector<
		boost::variant<double, std::shared_ptr<Matrix<double>>>
	>&data = vm->get_data();
	if (data.back().type() != typeid(double))return;
	double y = boost::get<double>(data.back());
	data.pop_back();
	if (data.back().type() != typeid(double))return;
	double &x = boost::get<double>(data.back());
	x = gamma_Q(x, y);
	return;
}
void factorial_builtin(client::vmachine  *vm) {
	std::vector<
		boost::variant<double, std::shared_ptr<Matrix<double>>>
	>&data = vm->get_data();
	if (data.back().type() != typeid(double))return;
	double &x = boost::get<double>(data.back());
	x = factorial(x);
	return;
}
void binomial_builtin(client::vmachine  *vm) {
	std::vector<
		boost::variant<double, std::shared_ptr<Matrix<double>>>
	>&data = vm->get_data();
	if (data.back().type() != typeid(double))return;
	double y = boost::get<double>(data.back());
	data.pop_back();
	if (data.back().type() != typeid(double))return;
	double &x = boost::get<double>(data.back());
	x = binomial(x, y);
	return;
}
void beta_builtin(client::vmachine  *vm) {
	std::vector<
		boost::variant<double, std::shared_ptr<Matrix<double>>>
	>&data = vm->get_data();
	if (data.back().type() != typeid(double))return;
	double y = boost::get<double>(data.back());
	data.pop_back();
	if (data.back().type() != typeid(double))return;
	double &x = boost::get<double>(data.back());
	x = beta(x, y);
	return;
}
void incomplete_beta_builtin(client::vmachine  *vm) {
	std::vector<
		boost::variant<double, std::shared_ptr<Matrix<double>>>
	>&data = vm->get_data();
	if (data.back().type() != typeid(double))return;
	double z = boost::get<double>(data.back());
	data.pop_back();
	if (data.back().type() != typeid(double))return;
	double y = boost::get<double>(data.back());
	data.pop_back();
	if (data.back().type() != typeid(double))return;
	double &x = boost::get<double>(data.back());
	x = incomplete_beta(x, y,z);
	return;
}
void error_function_builtin(client::vmachine  *vm) {
	std::vector<
		boost::variant<double, std::shared_ptr<Matrix<double>>>
	>&data = vm->get_data();
	if (data.back().type() != typeid(double))return;
	double &x = boost::get<double>(data.back());
	x = error_function(x);
	return;
}
void Legendre_p_builtin(client::vmachine  *vm) {
	std::vector<
		boost::variant<double, std::shared_ptr<Matrix<double>>>
	>&data = vm->get_data();
	if (data.back().type() != typeid(double))return;
	double y = boost::get<double>(data.back());
	data.pop_back();
	if (data.back().type() != typeid(double))return;
	double &x = boost::get<double>(data.back());
	x = Legendre_p(x, y);
	return;
}
void Legendre_q_builtin(client::vmachine  *vm) {
	std::vector<
		boost::variant<double, std::shared_ptr<Matrix<double>>>
	>&data = vm->get_data();
	if (data.back().type() != typeid(double))return;
	double y = boost::get<double>(data.back());
	data.pop_back();
	if (data.back().type() != typeid(double))return;
	double &x = boost::get<double>(data.back());
	x = Legendre_q(x, y);
	return;
}