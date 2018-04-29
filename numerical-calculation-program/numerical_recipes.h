#pragma once
#include<vector>
#include<string>
#include<memory>
#include"Matrix.h"
#include"Linear_solution.h"
#define EPS pow(2.0, -40.0)
using namespace std;

void Gauss_Jordan(Matrix<double> &A, Matrix<double> &b);//���A��������Ax=b
void solve_linear_LU_method(const Matrix<double> &A, Matrix<double> &b);//���Ax=b��AΪ����
std::vector<std::shared_ptr<Linear_solution<double>>> solve_linear_SVD_method(const Matrix<double> &A,const Matrix<double> &b, double eps = EPS);//Ax=bͨ�ýⷨ