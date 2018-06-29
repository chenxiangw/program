#include<iostream>
#include <cmath>
#include<limits>
#include"interpolation.h"
#include"Evaluation_function.hpp"
#define EPS 1.0e-10
#define ITER_MAX 20//���λ�����������
#define MID_ITER_MAX 15//�е㷨���λ�����������
using namespace std;

double trapezoidal_integration(double(*func)(double),const double a,const double b,const int n){
	/*
		���λ���
		������
			func����������
			a�����������½�
			b�����������Ͻ�
			n������2^(n-1)���ڵ����
		���أ�
			���ֽ��
	*/
	static double integration_value;
	if (n == 0) {
		return (integration_value = 0.5*(b - a)*((*func)(a) + (*func)(b)));//��ʼ������ֵ
	}
	else {
		int x_num = 1;
		for (int  j = 0; j<n - 1; ++j) x_num *=2;
		double interval = (b - a) / x_num;
		double x = a + 0.5*interval;
		double sum = 0.0;
		for (int j = 1; j <= x_num; ++j, x += interval) sum += (*func)(x);
		integration_value = 0.5*(integration_value + interval*sum);
		return integration_value;
	}
}

double trapezoidal_integration_iter(double(*func)(double), const double a, const double b){
	/*
	�������λ���
	������
	func����������
	a�����������½�
	b�����������Ͻ�
	���أ�
	���ֽ��
	*/
	double integration_value;
	double old = DBL_MAX;
	for (int j = 0; j < ITER_MAX; j++) {
		integration_value = trapezoidal_integration(func, a, b, j);
		if (fabs(integration_value - old) < EPS*fabs(old) && j > 4) return integration_value;
		old = integration_value;
	}
	cout << "���ֳ��������������ֵδ����" << endl;
	return 0.0;
}

double trapezoidal_integration_iter_Simpson(double(*func)(double), const double a, const double b){
	/*
	Simpson�������λ���
	������
	func����������
	a�����������½�
	b�����������Ͻ�
	���أ�
	���ֽ��
	*/
	double integration_value_weighted, integration_value;
	double old_weighted = DBL_MAX;
	double old = DBL_MAX;
	for (int j = 0; j < ITER_MAX; ++j) {
		integration_value = trapezoidal_integration(func, a, b, j);
		integration_value_weighted = (4.0*integration_value - old) / 3.0;
		if (fabs(integration_value_weighted - old_weighted) < EPS*fabs(old_weighted) && j > 4) return integration_value_weighted;
		old_weighted = integration_value_weighted;
		old = integration_value;
	}
	cout << "���ֳ��������������ֵδ����" << endl;
	return 0.0;
}

double trapezoidal_integration_iter_Romberg(double(*func)(double), const double a, const double b){
	/*
	Romberg�������λ���
	������
	func����������
	a�����������½�
	b�����������Ͻ�
	���أ�
	���ֽ��
	*/
	double integration_value;
	int order = 4;
	vector<double> s(ITER_MAX + 1), h(ITER_MAX + 1);
	h[0] = 1.0;
	for (int j = 0; j < ITER_MAX; ++j) {
		s[j] = trapezoidal_integration(func, a, b, j);
		if (j >= order) {
			double dy=polynomial_interpolation(h, s, 0.0, integration_value, j - order, j);
			if (fabs(dy) < EPS*fabs(integration_value)) return integration_value;
		}
		s[j + 1] = s[j];
		h[j + 1] = 0.25*h[j];
	}
	cout << "���ֳ��������������ֵδ����" << endl;
	return 0.0;
}

double mid_trapezoidal_integration(double(*func)(double), const double a, const double b, const int n) {
	/*
	�е㷨���λ��֣�������ֶ˵�ֵ
	������
	func����������
	a�����������½�
	b�����������Ͻ�
	n������3^(n-1)���ڵ����
	���أ�
	���ֽ��
	*/

	static double integration_value;
	if (n == 0) {
		return (integration_value = (b - a)*func(0.5*(a + b)));
	}
	else {
		int x_num = 1;
		for (int j = 0; j<n - 1; ++j) x_num *= 3;
		double interval = (b - a) / (3.0*x_num);
		double interval2 = interval * 2;
		double x = a + 0.5*interval;
		double sum = 0.0;
		for (int j = 1; j <= x_num; ++j) {
			sum += func(x);
			x += interval2;
			sum += func(x);
			x += interval;
		}
		integration_value = (integration_value + (b - a)*sum / x_num) / 3.0;
		return integration_value;
	}
}

double func_transformed_inf(double(*func)(double), double x) {
	return func(1.0 / x) / (x*x);
}

double mid_trapezoidal_integration_inf(double(*func)(double), const double a, const double b, const int n) {
	/*
	�е㷨���λ��֣�������ֶ˵�ֵ
	���ã�
		b�����a>0
		a�����b<0
	������
	func����������
	a�����������½�
	b�����������Ͻ�
	n������3^(n-1)���ڵ����
	���أ�
	���ֽ��
	*/

	double a_transformed = 1 / b;
	double b_transformed = 1 / a;
	static double integration_value;
	if (n == 0) {
		return (integration_value = (b_transformed - a_transformed)*func_transformed_inf(func, 0.5*(a_transformed + b_transformed)));
	}
	else {
		int x_num = 1;
		for (int j = 0; j<n - 1; ++j) x_num *= 3;
		double interval = (b_transformed - a_transformed) / (3.0*x_num);
		double interval2 = interval * 2;
		double x = a_transformed + 0.5*interval;
		double sum = 0.0;
		for (int j = 1; j <= x_num; ++j) {
			sum += func_transformed_inf(func,x);
			x += interval2;
			sum += func_transformed_inf(func, x);
			x += interval;
		}
		integration_value = (integration_value + (b_transformed - a_transformed)*sum / x_num) / 3.0;
		return integration_value;
	}
}

double func_transformed_singular_a(double(*func)(double),double a, double x) {
	return 2.0*x*func(a + x*x);
}

double mid_trapezoidal_integration_singular_a(double(*func)(double), const double a, const double b, const int n) {
	/*
	�е㷨���λ��֣�������ֶ˵�ֵ
	���ã�
		a��Ϊ���������
	������
		func����������
		a�����������½�
		b�����������Ͻ�
		n������3^(n-1)���ڵ����
	���أ�
		���ֽ��
	*/

	double a_transformed = 0.0;
	double b_transformed = sqrt(b-a);
	static double integration_value;
	if (n == 0) {
		return (integration_value = (b_transformed - a_transformed)*func_transformed_singular_a(func, a, 0.5*(a_transformed + b_transformed)));
	}
	else {
		int x_num = 1;
		for (int j = 0; j<n - 1; ++j) x_num *= 3;
		double interval = (b_transformed - a_transformed) / (3.0*x_num);
		double interval2 = interval * 2;
		double x = a_transformed + 0.5*interval;
		double sum = 0.0;
		for (int j = 1; j <= x_num; ++j) {
			sum += func_transformed_singular_a(func, a, x);
			x += interval2;
			sum += func_transformed_singular_a(func, a, x);
			x += interval;
		}
		integration_value = (integration_value + (b_transformed - a_transformed)*sum / x_num) / 3.0;
		return integration_value;
	}
}

double func_transformed_singular_b(double(*func)(double), double b, double x) {
	return 2.0*x*func(b - x*x);
}

double mid_trapezoidal_integration_singular_b(double(*func)(double), const double a, const double b, const int n) {
	/*
	�е㷨���λ��֣�������ֶ˵�ֵ
	���ã�
		b��Ϊ���������
	������
	func����������
	a�����������½�
	b�����������Ͻ�
	n������3^(n-1)���ڵ����
	���أ�
	���ֽ��
	*/

	double a_transformed = 0.0;
	double b_transformed = sqrt(b - a);
	static double integration_value;
	if (n == 0) {
		return (integration_value = (b_transformed - a_transformed)*func_transformed_singular_b(func, b, 0.5*(a_transformed + b_transformed)));
	}
	else {
		int x_num = 1;
		for (int j = 0; j<n - 1; ++j) x_num *= 3;
		double interval = (b_transformed - a_transformed) / (3.0*x_num);
		double interval2 = interval * 2;
		double x = a_transformed + 0.5*interval;
		double sum = 0.0;
		for (int j = 1; j <= x_num; ++j) {
			sum += func_transformed_singular_b(func, b, x);
			x += interval2;
			sum += func_transformed_singular_b(func, b, x);
			x += interval;
		}
		integration_value = (integration_value + (b_transformed - a_transformed)*sum / x_num) / 3.0;
		return integration_value;
	}
}

double mid_trapezoidal_integration_iter_Romberg(double(*func)(double), const double a, const double b,const int type=0){
	/*
	����������Ĺ������λ���
	������
	func����������
	a�����������½�
	b�����������Ͻ�
	type��
		0. ����������
		1. aΪ�����Ļ���
		2. bΪ�����Ļ���
	���أ�
	���ֽ��
	*/
	if (a >= b) {
		cout << "���ַ�Χ����" << endl;
		return 0.0;
	}
	double (*mid_trapezoidal_integration_chosen)(double(*)(double), const double, const double, const int)=nullptr;
	if (type == 0) {
		if(a == -DBL_MAX || b == DBL_MAX) mid_trapezoidal_integration_chosen = mid_trapezoidal_integration_inf;//����������
		else mid_trapezoidal_integration_chosen = mid_trapezoidal_integration;//��ͨ����
	}
	else if (type == 1) {
		mid_trapezoidal_integration_chosen = mid_trapezoidal_integration_singular_a;
	}
	else if (type == 2) {
		mid_trapezoidal_integration_chosen = mid_trapezoidal_integration_singular_b;
	}
	else {
		cout << "�������ʹ���" << endl;
		return 0.0;
	}
	double integration_value;
	int order = 4;
	vector<double> h(MID_ITER_MAX + 1), s(MID_ITER_MAX + 1);
	h[0] = 1.0;
	for (int j = 0; j < MID_ITER_MAX; ++j) {
		s[j] = mid_trapezoidal_integration_chosen(func, a, b, j);
		if (j >= order) {
			double dy = polynomial_interpolation(h, s, 0.0, integration_value, j - order, j);
			if (fabs(dy) < EPS*fabs(integration_value)) return integration_value;
		}
		s[j + 1] = s[j];
		h[j + 1] = h[j] / 9.0;
	}
	cout << "���ֳ��������������ֵδ����" << endl;
	return 0.0;
}

void Chebyshev_integration(const double a, const double b, const vector<double> &c, vector<double> &c_integrate) {
	/*
	�б�ѩ�����ʽ�ĵ���
	������
	a�������½�
	b�������Ͻ�
	c���б�ѩ�����ʽϵ��
	c_integrate���б�ѩ�����ʽ���ֵ�ϵ��
	*/
	int n = c.size();
	double sum = 0.0, sign = 1.0;
	double factor = 0.25*(b - a);
	for (int j = 1; j <= n - 2; j++) {
		c_integrate[j] = factor * (c[j - 1] - c[j + 1]) / j;
		sum += sign * c_integrate[j];
		sign = -sign;
	}
	c_integrate[n - 1] = factor * c[n - 2] / (n - 1);
	sum += sign * c_integrate[n - 1];
	c_integrate[0] = 2.0*sum;
}

double mid_trapezoidal_integration(Evaluation_function &func, const double a, const double b, const int n) {
	/*
	�е㷨���λ��֣�������ֶ˵�ֵ
	������
	func����������
	a�����������½�
	b�����������Ͻ�
	n������3^(n-1)���ڵ����
	���أ�
	���ֽ��
	*/

	static double integration_value;
	if (n == 0) {
		return (integration_value = (b - a)*func(0.5*(a + b)));
	}
	else {
		int x_num = 1;
		for (int j = 0; j<n - 1; ++j) x_num *= 3;
		double interval = (b - a) / (3.0*x_num);
		double interval2 = interval * 2;
		double x = a + 0.5*interval;
		double sum = 0.0;
		for (int j = 1; j <= x_num; ++j) {
			sum += func(x);
			x += interval2;
			sum += func(x);
			x += interval;
		}
		integration_value = (integration_value + (b - a)*sum / x_num) / 3.0;
		return integration_value;
	}
}

double func_transformed_inf(Evaluation_function &func, double x) {
	return func(1.0 / x) / (x*x);
}

double mid_trapezoidal_integration_inf(Evaluation_function &func, const double a, const double b, const int n) {
	/*
	�е㷨���λ��֣�������ֶ˵�ֵ
	���ã�
	b�����a>0
	a�����b<0
	������
	func����������
	a�����������½�
	b�����������Ͻ�
	n������3^(n-1)���ڵ����
	���أ�
	���ֽ��
	*/

	double a_transformed = 1 / b;
	double b_transformed = 1 / a;
	static double integration_value;
	if (n == 0) {
		return (integration_value = (b_transformed - a_transformed)*func_transformed_inf(func, 0.5*(a_transformed + b_transformed)));
	}
	else {
		int x_num = 1;
		for (int j = 0; j<n - 1; ++j) x_num *= 3;
		double interval = (b_transformed - a_transformed) / (3.0*x_num);
		double interval2 = interval * 2;
		double x = a_transformed + 0.5*interval;
		double sum = 0.0;
		for (int j = 1; j <= x_num; ++j) {
			sum += func_transformed_inf(func, x);
			x += interval2;
			sum += func_transformed_inf(func, x);
			x += interval;
		}
		integration_value = (integration_value + (b_transformed - a_transformed)*sum / x_num) / 3.0;
		return integration_value;
	}
}

double func_transformed_singular_a(Evaluation_function &func, double a, double x) {
	return 2.0*x*func(a + x*x);
}

double mid_trapezoidal_integration_singular_a(Evaluation_function &func, const double a, const double b, const int n) {
	/*
	�е㷨���λ��֣�������ֶ˵�ֵ
	���ã�
	a��Ϊ���������
	������
	func����������
	a�����������½�
	b�����������Ͻ�
	n������3^(n-1)���ڵ����
	���أ�
	���ֽ��
	*/

	double a_transformed = 0.0;
	double b_transformed = sqrt(b - a);
	static double integration_value;
	if (n == 0) {
		return (integration_value = (b_transformed - a_transformed)*func_transformed_singular_a(func, a, 0.5*(a_transformed + b_transformed)));
	}
	else {
		int x_num = 1;
		for (int j = 0; j<n - 1; ++j) x_num *= 3;
		double interval = (b_transformed - a_transformed) / (3.0*x_num);
		double interval2 = interval * 2;
		double x = a_transformed + 0.5*interval;
		double sum = 0.0;
		for (int j = 1; j <= x_num; ++j) {
			sum += func_transformed_singular_a(func, a, x);
			x += interval2;
			sum += func_transformed_singular_a(func, a, x);
			x += interval;
		}
		integration_value = (integration_value + (b_transformed - a_transformed)*sum / x_num) / 3.0;
		return integration_value;
	}
}

double func_transformed_singular_b(Evaluation_function &func, double b, double x) {
	return 2.0*x*func(b - x*x);
}

double mid_trapezoidal_integration_singular_b(Evaluation_function &func, const double a, const double b, const int n) {
	/*
	�е㷨���λ��֣�������ֶ˵�ֵ
	���ã�
	b��Ϊ���������
	������
	func����������
	a�����������½�
	b�����������Ͻ�
	n������3^(n-1)���ڵ����
	���أ�
	���ֽ��
	*/

	double a_transformed = 0.0;
	double b_transformed = sqrt(b - a);
	static double integration_value;
	if (n == 0) {
		return (integration_value = (b_transformed - a_transformed)*func_transformed_singular_b(func, b, 0.5*(a_transformed + b_transformed)));
	}
	else {
		int x_num = 1;
		for (int j = 0; j<n - 1; ++j) x_num *= 3;
		double interval = (b_transformed - a_transformed) / (3.0*x_num);
		double interval2 = interval * 2;
		double x = a_transformed + 0.5*interval;
		double sum = 0.0;
		for (int j = 1; j <= x_num; ++j) {
			sum += func_transformed_singular_b(func, b, x);
			x += interval2;
			sum += func_transformed_singular_b(func, b, x);
			x += interval;
		}
		integration_value = (integration_value + (b_transformed - a_transformed)*sum / x_num) / 3.0;
		return integration_value;
	}
}

double mid_trapezoidal_integration_iter_Romberg(Evaluation_function &func, const double a, const double b, const int type = 0) {
	/*
	����������Ĺ������λ���
	������
	func����������
	a�����������½�
	b�����������Ͻ�
	type��
	0. ����������
	1. aΪ�����Ļ���
	2. bΪ�����Ļ���
	���أ�
	���ֽ��
	*/
	if (a >= b) {
		cout << "���ַ�Χ����" << endl;
		return 0.0;
	}
	double(*mid_trapezoidal_integration_chosen)(Evaluation_function&, const double, const double, const int) = nullptr;
	if (type == 0) {
		if (a == -DBL_MAX || b == DBL_MAX) mid_trapezoidal_integration_chosen = mid_trapezoidal_integration_inf;//����������
		else mid_trapezoidal_integration_chosen = mid_trapezoidal_integration;//��ͨ����
	}
	else if (type == 1) {
		mid_trapezoidal_integration_chosen = mid_trapezoidal_integration_singular_a;
	}
	else if (type == 2) {
		mid_trapezoidal_integration_chosen = mid_trapezoidal_integration_singular_b;
	}
	else {
		cout << "�������ʹ���" << endl;
		return 0.0;
	}
	double integration_value;
	int order = 4;
	vector<double> h(MID_ITER_MAX + 1), s(MID_ITER_MAX + 1);
	h[0] = 1.0;
	for (int j = 0; j < MID_ITER_MAX; ++j) {
		s[j] = mid_trapezoidal_integration_chosen(func, a, b, j);
		if (j >= order) {
			double dy = polynomial_interpolation(h, s, 0.0, integration_value, j - order, j);
			if (fabs(dy) < EPS*fabs(integration_value)) return integration_value;
		}
		s[j + 1] = s[j];
		h[j + 1] = h[j] / 9.0;
	}
	cout << "���ֳ��������������ֵδ����" << endl;
	return 0.0;
}



#undef EPS
#undef ITER_MAX
#undef MID_ITER_MAX