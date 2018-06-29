#pragma once
#include<iostream>
#include<string>
#include<map>
#include<utility>
#include<vector>
#include<limits>
//#include"Evaluation_function.h"
#include"core/numerical_recipes.h"

#include"parser/vm.hpp"

void inverse_builtin(client::vmachine  *vm);
void determinant_builtin(client::vmachine  *vm);
void solve_builtin(client::vmachine  *vm);

//void integrate(const std::vector<std::string> *input);
//void diff(const std::vector<std::string> *input);

void cos_builtin(client::vmachine  *vm);
void sin_builtin(client::vmachine  *vm);
void tan_builtin(client::vmachine  *vm);
void acos_builtin(client::vmachine  *vm);
void asin_builtin(client::vmachine  *vm);
void atan_builtin(client::vmachine  *vm);

void gamma_builtin(client::vmachine  *vm);
void ln_gamma_builtin(client::vmachine  *vm);
void gamma_P_builtin(client::vmachine  *vm);
void gamma_Q_builtin(client::vmachine  *vm);
void factorial_builtin(client::vmachine  *vm);
void binomial_builtin(client::vmachine  *vm);
void beta_builtin(client::vmachine  *vm);
void incomplete_beta_builtin(client::vmachine  *vm);
void error_function_builtin(client::vmachine  *vm);
void Legendre_p_builtin(client::vmachine  *vm);
void Legendre_q_builtin(client::vmachine  *vm);