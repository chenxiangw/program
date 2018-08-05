#pragma once
#include<vector>
#include<string>
#include<memory>
#include"Matrix.hpp"
#include"Linear_solution.hpp"
#include"Evaluation_function.hpp"
#define EPS pow(2.0, -100.0)


void Gauss_Jordan(Matrix<double> &A, Matrix<double> &b);//���A��������Ax=b
void solve_linear_LU_method(const Matrix<double> &A, Matrix<double> &b);//���Ax=b��AΪ����
std::vector<std::shared_ptr<Linear_solution<double>>> solve_linear_SVD_method(const Matrix<double> &A,const Matrix<double> &b, double eps = EPS);//Ax=bͨ�ýⷨ
double mid_trapezoidal_integration_iter_Romberg(Evaluation_function &func, const double a, const double b, const int type = 0);
double derivation_Ridders(Evaluation_function &func, double &y1order, const double x, const double init_step);
