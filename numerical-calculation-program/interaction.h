#pragma once
#include<string>
#include<unordered_map>
#include <functional> 
#include"Matrix.h"

struct Function {};

template<typename ...Arguments>
struct Function_type :public Function {
	std::function<void(Arguments...)> callback;
	Function_type(std::function<void(Arguments...)> input) : callback(input) {}
};

void interaction();
void dict(const std::string &input);
void matrix(const std::string &input);
void clear(const std::string &input);
void inverse(const std::string &input);
void determinant(const std::string &input);
void solve(const std::string &input);