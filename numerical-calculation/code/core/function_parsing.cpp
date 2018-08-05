#include<string>
#include<iostream>
#include<cmath>
#include<vector>
#include <boost\spirit\include\qi.hpp>
#include<boost\spirit\include\phoenix.hpp>
#include<boost\bind\bind.hpp>
#include <boost\lambda\lambda.hpp>

void operator_cos(std::vector<double> *nums) {
	int n = nums->size();
	if (n < 1) {
		throw  std::string("表达式错误");
	}
	(*nums)[n - 1] = cos((*nums)[n - 1]);
}
void operator_sin(std::vector<double> *nums) {
	int n = nums->size();
	if (n < 1) {
		throw  std::string("表达式错误");
	}
	(*nums)[n - 1] = sin((*nums)[n - 1]);
}
void operator_tan(std::vector<double> *nums) {
	int n = nums->size();
	if (n < 1) {
		throw  std::string("表达式错误");
	}
	(*nums)[n - 1] = tan((*nums)[n - 1]);
}
void operator_acos(std::vector<double> *nums) {
	int n = nums->size();
	if (n < 1) {
		throw  std::string("表达式错误");
	}
	(*nums)[n - 1] = acos((*nums)[n - 1]);
}
void operator_asin(std::vector<double> *nums) {
	int n = nums->size();
	if (n < 1) {
		throw  std::string("表达式错误");
	}
	(*nums)[n - 1] = asin((*nums)[n - 1]);
}
void operator_atan(std::vector<double> *nums) {
	int n = nums->size();
	if (n < 1) {
		throw  std::string("表达式错误");
	}
	(*nums)[n - 1] = atan((*nums)[n - 1]);
}

void operator_pow(std::vector<double> *nums) {
	int n = nums->size();
	if (n < 2) {
		throw  std::string("表达式错误");
	}
	(*nums)[n - 2] = pow((*nums)[n - 2], (*nums)[n - 1]);
	nums->pop_back();
}
void operator_mutiply(std::vector<double> *nums) {
	int n = nums->size();
	if (n < 2) {
		throw  std::string("表达式错误");
	}
	(*nums)[n - 2] *= (*nums)[n - 1];
	nums->pop_back();
}
void operator_divide(std::vector<double> *nums) {
	int n = nums->size();
	if (n < 2) {
		throw  std::string("表达式错误");
	}
	if ((*nums)[n - 1] == 0) {
		throw  std::string("表达式除以0");
	}
	(*nums)[n - 2] /= (*nums)[n - 1];
	nums->pop_back();
}
void operator_add(std::vector<double> *nums) {
	int n = nums->size();
	if (n < 2) {
		throw  std::string("表达式错误");
	}
	(*nums)[n - 2] += (*nums)[n - 1];
	nums->pop_back();
}
void operator_minus(std::vector<double> *nums) {
	int n = nums->size();
	if (n < 2) {
		throw  std::string("表达式错误");
	}
	(*nums)[n - 2] -= (*nums)[n - 1];
	nums->pop_back();
}



bool function_parsing(std::string &function, const double value, double &result, const std::string var = "x") {
	using namespace boost::spirit::qi;
	using boost::phoenix::ref;
	using boost::phoenix::push_back;
	using boost::bind;
	auto begin = function.begin();
	auto end = function.end();
	rule<std::string::iterator, space_type, double()> priority0,//函数运算、读取数值和括号操作
		priority1,//幂运算
		priority2,//乘除法运算
		priority3;//加减法运算

	std::vector<double> nums;
	priority0 = 
		string("cos") >> '(' >> priority3[bind(operator_cos, &nums)] >> ')' |
		string("sin") >> '(' >> priority3[bind(operator_sin, &nums)] >> ')' |
		string("tan") >> '(' >> priority3[bind(operator_tan, &nums)] >> ')' |
		string("acos") >> '(' >> priority3[bind(operator_acos, &nums)] >> ')' |
		string("asin") >> '(' >> priority3[bind(operator_asin, &nums)] >> ')' |
		string("atan") >> '(' >> priority3[bind(operator_atan, &nums)] >> ')' |
		string(var)[push_back(boost::phoenix::ref(nums), value)] | double_[push_back(boost::phoenix::ref(nums), _1)] | ('(' >> priority3 >> ')');

	priority1 = priority0 >> *('^' >> priority0[bind(operator_pow, &nums)]);
	priority2 = priority1 >> *(('*' >> priority1[bind(operator_mutiply, &nums)]) | ('/' >> priority1[bind(operator_divide, &nums)]));
	priority3 = priority2 >> *(('+' >> priority2[bind(operator_add, &nums)]) | ('-' >> priority2[bind(operator_minus, &nums)]));
	bool succeed = phrase_parse(begin, end, priority3, space);
	if (begin != end)
		return false;
	if (nums.size() == 1) {
		result = nums[0];
	}
	else {
		std::cout << "函数错误" << std::endl;
		return false;
	}
	return succeed;
}
