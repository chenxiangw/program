#pragma once
#include<iostream>
#include<string>
#include<map>
#include<utility>
#include<vector>
#include<limits>
#include <boost/variant/recursive_variant.hpp>
#include"core/numerical_recipes.h"
#include"parser/program.hpp"


bool inverse_builtin(std::vector<data_type>&data);
bool determinant_builtin(std::vector<data_type>&data);
bool solve_builtin(std::vector<data_type>&data);

//bool integrate(const std::vector<std::string> *input);
//bool diff(const std::vector<std::string> *input);

bool cos_builtin(std::vector<data_type>&data);
bool sin_builtin(std::vector<data_type>&data);
bool tan_builtin(std::vector<data_type>&data);
bool acos_builtin(std::vector<data_type>&data);
bool asin_builtin(std::vector<data_type>&data);
bool atan_builtin(std::vector<data_type>&data);

bool gamma_builtin(std::vector<data_type>&data);
bool ln_gamma_builtin(std::vector<data_type>&data);
bool gamma_P_builtin(std::vector<data_type>&data);
bool gamma_Q_builtin(std::vector<data_type>&data);
bool factorial_builtin(std::vector<data_type>&data);
bool binomial_builtin(std::vector<data_type>&data);
bool beta_builtin(std::vector<data_type>&data);
bool incomplete_beta_builtin(std::vector<data_type>&data);
bool error_function_builtin(std::vector<data_type>&data);
bool Legendre_p_builtin(std::vector<data_type>&data);
bool Legendre_q_builtin(std::vector<data_type>&data);