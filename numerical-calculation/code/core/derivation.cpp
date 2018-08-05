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
	Ridders����ʽ���Ʒ���
	������
	func����ֵ����
	x����λ��
	init_step����ʼ����
	y1order����õ�1�׵���
	���أ�
	���
	*/
	const double step_decrease = 1.4, step_decrease2 = (step_decrease*step_decrease);
	const int NTAB = 10;//���ƽ���
	const double SAFE = 2.0;//��ȫϵ��

	if (fabs(init_step) < EPS) {
		cout << "�󵼴��󣺲�������Ϊ0" << endl;
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
	�б�ѩ�����ʽ�ĵ���
	������
	a�������½�
	b�������Ͻ�
	c���б�ѩ�����ʽϵ��
	c_diff���б�ѩ�����ʽ���ֵ�ϵ��
	*/

	size_t n = c.size();
	c_diff[n - 1] = 0.0;
	c_diff[n - 2] = 2 * (n - 1)*c[n - 1];
	for (size_t j = n - 3; j >= 0; j--)c_diff[j] = c_diff[j + 2] + 2 * (j + 1)*c[j + 1];
	double factor = 2.0 / (b - a);
	for (int j = 0; j<n; j++)
		c_diff[j] *= factor;
}
