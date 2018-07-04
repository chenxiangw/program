

#include<cmath>

#include <boost/spirit/include/qi.hpp>

#include "builtin_functions.h"
#include "core/special_functions.h"
#include "core/numerical_recipes.h"


bool inverse_builtin(std::vector<data_type>&data) {
	if (data.size() != 1)return false;
	if (data.back().type() != typeid(std::shared_ptr<Matrix<double>>))return false;
	std::shared_ptr<Matrix<double>> &matrix = boost::get< std::shared_ptr<Matrix<double>>>(data.back());
	try {
		auto inverse = matrix->inverse();
		matrix = inverse;
	}
	catch (std::string &error) {
		std::cout << error << std::endl;;
	}
	return true;
}
bool determinant_builtin(std::vector<data_type>&data){
	if (data.size() != 1)return false;
	if (data.back().type() != typeid(std::shared_ptr<Matrix<double>>))return false;
	std::shared_ptr<Matrix<double>> &matrix = boost::get< std::shared_ptr<Matrix<double>>>(data.back());
	try {
		std::shared_ptr<double> result = matrix->determinant();
		data.pop_back();
		if (result != nullptr)data.push_back(result);
	}
	catch (std::string &error) {
		std::cout << error << std::endl;;
	}
	return true;
}
bool solve_builtin(std::vector<data_type>&data) {
	if (data.size() != 2)return false;
	if (data.back().type() != typeid(std::shared_ptr<Matrix<double>>))return false;
	std::shared_ptr<Matrix<double>> b = boost::get< std::shared_ptr<Matrix<double>>>(data.back());
	data.pop_back();
	std::shared_ptr<Matrix<double>> A = boost::get< std::shared_ptr<Matrix<double>>>(data.back());
	data.pop_back();
	if (A->msize() != b->msize()) {
		std::cout << "矩阵格式错误" << std::endl;
		return false;
	}
	auto result = solve_linear_SVD_method(*A, *b);
	size_t question_num = result.size();
	for (int i = 0; i < question_num; ++i) {
		std::cout << "第" << i + 1 << "个方程：" << std::endl;
		result[i]->print();
	}
	return true;
}

//bool integrate(const std::vector<std::string> *input) {
//	using namespace std;
//	if (input->size() != 3) {
//		cout << "参数错误" << endl;
//		return false;
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
//	return false;
//}
//bool diff(const std::vector<std::string> *input) {
//	using namespace std;
//	if (input->size() != 2) {
//		cout << "参数错误" << endl;
//		return false;
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
//	return false;
//}

bool cos_builtin(std::vector<data_type>&data) {
	if (data.size() != 1)return false;
	if (data.back().type() != typeid(std::shared_ptr<double>))return false;
	std::shared_ptr<double> &x = boost::get<std::shared_ptr<double>>(data.back());
	*x = cos(*x);
	return true;
}
bool sin_builtin(std::vector<data_type>&data) {
	if (data.size() != 1)return false;
	if (data.back().type() != typeid(std::shared_ptr<double>))return false;
	std::shared_ptr<double> &x = boost::get<std::shared_ptr<double>>(data.back());
	*x = sin(*x);
	return true;
}
bool tan_builtin(std::vector<data_type>&data) {
	if (data.size() != 1)return false;
	if (data.back().type() != typeid(std::shared_ptr<double>))return false;
	std::shared_ptr<double> &x = boost::get<std::shared_ptr<double>>(data.back());
	*x = tan(*x);
	return true;
}
bool acos_builtin(std::vector<data_type>&data) {
	if (data.size() != 1)return false;
	if (data.back().type() != typeid(std::shared_ptr<double>))return false;
	std::shared_ptr<double> &x = boost::get<std::shared_ptr<double>>(data.back());
	*x = acos(*x);
	return true;
}
bool asin_builtin(std::vector<data_type>&data) {
	if (data.size() != 1)return false;
	if (data.back().type() != typeid(std::shared_ptr<double>))return false;
	std::shared_ptr<double> &x = boost::get<std::shared_ptr<double>>(data.back());
	*x = asin(*x);
	return true;
}
bool atan_builtin(std::vector<data_type>&data) {
	if (data.size() != 1)return false;
	if (data.back().type() != typeid(std::shared_ptr<double>))return false;
	std::shared_ptr<double> &x = boost::get<std::shared_ptr<double>>(data.back());
	*x = atan(*x);
	return true;
}

bool gamma_builtin(std::vector<data_type>&data) {
	if (data.size() != 1)return false;
	if (data.back().type() != typeid(std::shared_ptr<double>))return false;
	std::shared_ptr<double> &x = boost::get<std::shared_ptr<double>>(data.back());
	*x = gamma(*x);
	return true;
}
bool ln_gamma_builtin(std::vector<data_type>&data) {
	if (data.size() != 1)return false;
	if (data.back().type() != typeid(std::shared_ptr<double>))return false;
	std::shared_ptr<double> &x = boost::get<std::shared_ptr<double>>(data.back());
	*x = ln_gamma(*x);
	return true;
}
bool gamma_P_builtin(std::vector<data_type>&data) {
	if (data.size() != 2)return false;
	if (data.back().type() != typeid(std::shared_ptr<double>))return false;
	std::shared_ptr<double> y = boost::get<std::shared_ptr<double>>(data.back());
	data.pop_back();
	if (data.back().type() != typeid(std::shared_ptr<double>))return false;
	std::shared_ptr<double> &x = boost::get<std::shared_ptr<double>>(data.back());
	*x = gamma_P(*x, *y);
	return true;
}
bool gamma_Q_builtin(std::vector<data_type>&data) {
	if (data.size() != 2)return false;
	if (data.back().type() != typeid(std::shared_ptr<double>))return false;
	std::shared_ptr<double> y = boost::get<std::shared_ptr<double>>(data.back());
	data.pop_back();
	if (data.back().type() != typeid(std::shared_ptr<double>))return false;
	std::shared_ptr<double> &x = boost::get<std::shared_ptr<double>>(data.back());
	*x = gamma_Q(*x, *y);
	return true;
}
bool factorial_builtin(std::vector<data_type>&data) {
	if (data.size() != 1)return false;
	if (data.back().type() != typeid(std::shared_ptr<double>))return false;
	std::shared_ptr<double> &x = boost::get<std::shared_ptr<double>>(data.back());
	*x = factorial(*x);
	return true;
}
bool binomial_builtin(std::vector<data_type>&data) {
	if (data.size() != 2)return false;
	if (data.back().type() != typeid(std::shared_ptr<double>))return false;
	std::shared_ptr<double> y = boost::get<std::shared_ptr<double>>(data.back());
	data.pop_back();
	if (data.back().type() != typeid(std::shared_ptr<double>))return false;
	std::shared_ptr<double> &x = boost::get<std::shared_ptr<double>>(data.back());
	*x = binomial(*x, *y);
	return true;
}
bool beta_builtin(std::vector<data_type>&data) {
	if (data.size() != 2)return false;
	if (data.back().type() != typeid(std::shared_ptr<double>))return false;
	std::shared_ptr<double> y = boost::get<std::shared_ptr<double>>(data.back());
	data.pop_back();
	if (data.back().type() != typeid(std::shared_ptr<double>))return false;
	std::shared_ptr<double> &x = boost::get<std::shared_ptr<double>>(data.back());
	*x = beta(*x, *y);
	return true;
}
bool incomplete_beta_builtin(std::vector<data_type>&data) {
	if (data.size() != 3)return false;
	if (data.back().type() != typeid(std::shared_ptr<double>))return false;
	std::shared_ptr<double> z = boost::get<std::shared_ptr<double>>(data.back());
	data.pop_back();
	if (data.back().type() != typeid(std::shared_ptr<double>))return false;
	std::shared_ptr<double> y = boost::get<std::shared_ptr<double>>(data.back());
	data.pop_back();
	if (data.back().type() != typeid(std::shared_ptr<double>))return false;
	std::shared_ptr<double> &x = boost::get<std::shared_ptr<double>>(data.back());
	*x = incomplete_beta(*x, *y, *z);
	return true;
}
bool error_function_builtin(std::vector<data_type>&data) {
	if (data.size() != 1)return false;
	if (data.back().type() != typeid(std::shared_ptr<double>))return false;
	std::shared_ptr<double> &x = boost::get<std::shared_ptr<double>>(data.back());
	*x = error_function(*x);
	return true;
}
bool Legendre_p_builtin(std::vector<data_type>&data) {
	if (data.size() != 2)return false;
	if (data.back().type() != typeid(std::shared_ptr<double>))return false;
	std::shared_ptr<double> y = boost::get<std::shared_ptr<double>>(data.back());
	data.pop_back();
	if (data.back().type() != typeid(std::shared_ptr<double>))return false;
	std::shared_ptr<double> &x = boost::get<std::shared_ptr<double>>(data.back());
	*x = Legendre_p(*x, *y);
	return true;
}
bool Legendre_q_builtin(std::vector<data_type>&data) {
	if (data.size() !=2)return false;
	if (data.back().type() != typeid(std::shared_ptr<double>))return false;
	std::shared_ptr<double> y = boost::get<std::shared_ptr<double>>(data.back());
	data.pop_back();
	if (data.back().type() != typeid(std::shared_ptr<double>))return false;
	std::shared_ptr<double> &x = boost::get<std::shared_ptr<double>>(data.back());
	*x = Legendre_q(*x, *y);
	return true;
}