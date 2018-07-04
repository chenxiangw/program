#include <cmath>
#include<iostream>
#include<limits>
#include<vector>
#include"Evaluation_function.hpp"
#include"numerical_recipes.h"
#include"Matrix.hpp"
#include"interpolation.h"
using namespace std;

double derivation_Ridders(Evaluation_function &func, double &y1order,const double x,const double init_step){
	/*
	Ridders多项式外推法求导
	参数：
	func：求值函数
	x：求导位置
	init_step：初始步长
	y1order：求得的1阶导数
	返回：
	误差
	*/
	const double step_decrease = 1.4, step_decrease2 = (step_decrease*step_decrease);
	const int NTAB = 10;//外推阶数
	const double SAFE = 2.0;//误差安全系数

	if (fabs(init_step) < EPS) {
		cout << "求导错误：步长不能为0" << endl;
		return 0;
	}
	Matrix<double> table = Matrix<double>(NTAB, NTAB);
	double step = init_step;
	table[0][0] = (func(x + step) - func(x - step))/ (2.0*step);
	double error = DBL_MAX;
	for (int i = 1; i < NTAB; ++i) {
		step /= step_decrease;
		table[0][i] = (func(x + step) - func(x - step)) / (2.0*step);
		double factor = step_decrease2;
		for (int j = 1; j <= i; ++j) {
			table[j][i] = (table[j - 1][i] * factor - table[j - 1][i - 1]) / (factor - 1.0);
			factor = step_decrease2*factor;
			double error_temp = max(fabs(table[j][i] - table[j - 1][i]), fabs(table[j][i] - table[j - 1][i - 1]));
			if (error_temp <= error) {
				error = error_temp;
				y1order = table[j][i];
			}
		}
		if (fabs(table[i][i] - table[i - 1][i - 1]) >= SAFE*error) break;
	}
	return error;
}

void Chebyshev_derivation(const double a, const double b, const vector<double> &c, vector<double> &c_diff){
	/*
	切比雪夫多项式的导数
	参数：
	a：区间下界
	b：区间上界
	c：切比雪夫多项式系数
	c_diff：切比雪夫多项式积分的系数
	*/

	size_t n = c.size();
	c_diff[n - 1] = 0.0;
	c_diff[n - 2] = 2 * (n - 1)*c[n - 1];
	for (size_t j = n - 3; j >= 0; j--)c_diff[j] = c_diff[j + 2] + 2 * (j + 1)*c[j + 1];
	double factor = 2.0 / (b - a);
	for (int j = 0; j<n; j++)
		c_diff[j] *= factor;
}
