#pragma once
#include<string>
#include<iostream>
#include<map>
extern std::map < std::string, std::string > functions;
bool function_parsing(std::string &function, const double value, double &result, const std::string var = "x");

class Evaluation_function {
	std::string function;
public:
	Evaluation_function(std::string function_name) {
		function = functions[function_name];
	}
	double operator()(double x) {
		double result = 0.0;
		bool succeed = function_parsing(function, x, result);
		if (!succeed) {
			throw  std::string("º¯ÊýÇóÖµ´íÎó");
		}
		return result;
	}
};
