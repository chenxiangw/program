#include<vector>
#include<iostream>
#include<memory>
#include"numerical_recipes.h"
using namespace std;

double polynomial_interpolation(const vector<double> &x_values,const vector<double> &y_values, const double x, double &y,const int start=0,const int end=-1) {
	/*
		����ʽ��ֵ����
		������
			x_values����֪xֵ
			y_values����֪��Ӧyֵ
			x����Ҫ�õ�yֵ��xλ��
			y����������õ���yֵ
			start��x��ʼ���λ��
			end��x�������λ��
		���أ�
			������dy
	*/
	int n;
	if(end==-1)n = x_values.size();
	else n = end - start;
	vector<double> C(n), D(n);
	//Ѱ��x_values����xֵ����ĵ���Ϊ��ʼ�ƽ�ֵ
	double diff = fabs(x - x_values[start]);
	int closet_x = start;
	for (int i = start; i < end; ++i) {
		double diff_temp = fabs(x - x_values[i]);
		if (diff_temp < 0) {
			y = y_values[i];
			return 0;
		}
		else if (diff_temp < diff) {
			closet_x = i;
			diff = diff_temp;
		}
		C[i-start] = y_values[i];
		D[i-start] = y_values[i];
	}
	y = y_values[closet_x--];
	double dy = 0;
	for (int m = 1; m < n; ++m) {
		for (int i = start; i < n - m+start; ++i) {
			double C_factor = x_values[i] - x;
			double D_factor = x_values[i+m] - x;
			double denominator = C_factor - D_factor;
			if (fabs(denominator) < EPS) {
				cout << "��ֵ�������������ͬx��" << endl;
				return 0;
			}
			denominator = (C[i + 1-start] - D[i-start]) / denominator;
			C[i-start] = C_factor*denominator;
			D[i-start] = D_factor*denominator;
		}
		dy = (2 * (closet_x + 1)) < (n - m) ? C[closet_x + 1-start] : D[(--closet_x)+1-start];
		y += dy;
	}
	return dy;
}

double rational_interpolation(const vector<double> &x_values, const vector<double> &y_values, const double x, double &y) {
	/*
		�����ֵ����
		������
		x_values����֪xֵ
		y_values����֪��Ӧyֵ
		x����Ҫ�õ�yֵ��xλ��
		y����������õ���yֵ
		���أ�
		������dy
	*/
	int n = x_values.size();
	vector<double> C(n), D(n);
	//Ѱ��x_values����xֵ����ĵ���Ϊ��ʼ�ƽ�ֵ
	double diff = fabs(x - x_values[0]);
	int closet_x = 0;
	for (int i = 0; i < n; ++i) {
		double diff_temp = fabs(x - x_values[i]);
		if (diff_temp < 0) {
			y = y_values[i];
			return 0;
		}
		else if (diff_temp < diff) {
			closet_x = i;
			diff = diff_temp;
		}
		C[i] = y_values[i];
		D[i] = y_values[i];
	}
	y = y_values[closet_x--];
	double dy = 0;
	for (int m = 1; m < n; ++m) {
		for (int i = 0; i < (n - m); ++i) {
			double C_factor = (x_values[i] - x)*D[i] / (x_values[i + m] - x);
			double both_factor = C_factor - C[i + 1];
			if (fabs(both_factor) < EPS) {
				cout << "��ֵ������ֵ��Ϊ����" << endl;
				return 0;
			}
			both_factor = (C[i + 1] - D[i]) / both_factor;
			D[i] = C[i + 1] * both_factor;
			C[i] = C_factor*both_factor;
		}
		dy = (2 * (closet_x + 1)) < (n - m) ? C[closet_x + 1] : D[closet_x--];
		y += dy;
	}
	return dy;
}

shared_ptr<vector<double>> spline_solve_y2order(const vector<double> &x_values, const vector<double> &y_values, const double y_boundary1,const double y_boundary2) {
	/*
		����������ֵ
		������
		x_values����֪xֵ
		y_values����֪��Ӧyֵ
		y_boundary1����һ���㴦��һ�׵���
		y_boundary2�����һ���㴦��һ�׵���
		���أ�
		��ֵ����ʽ�Ķ��׵���
	*/
	int n = x_values.size();
	vector<double> temp(n);
	shared_ptr<vector<double>> y2order(new vector<double>(n));
	if (fabs(y_boundary1) > pow(2, 100)) (*y2order)[0] = temp[0] = 0.0;//�±߽絼��ֵ�����������Ȼ�߽�����
	else {
		(*y2order)[0] = -0.5;
		temp[0] = (3.0 / (x_values[1] - x_values[0]))*((y_values[1] - y_values[0]) / (x_values[1] - x_values[0]) - y_boundary1);
	}
	for (int i = 1; i < n - 1; ++i) {
		double sig = (x_values[i] - x_values[i - 1]) / (x_values[i + 1] - x_values[i - 1]);
		double p = sig*(*y2order)[i - 1] + 2.0;
		(*y2order)[i] = (sig - 1.0) / p;
		temp[i] = (y_values[i + 1] - y_values[i]) / (x_values[i + 1] - x_values[i]) - (y_values[i] - y_values[i - 1]) / (x_values[i] - x_values[i - 1]);
		temp[i] = (6.0*temp[i] / (x_values[i + 1] - x_values[i - 1]) - sig*temp[i - 1]) / p;
	}
	if (fabs(y_boundary2) > pow(2, 100)) (*y2order)[n - 1] = temp[n - 1] = 0.0;
	else {
		(*y2order)[n - 1] = 0.5;
		temp[n-1]= (3.0 / (x_values[n-1] - x_values[n-2]))*(y_boundary2-(y_values[n-1] - y_values[n-2]) / (x_values[n-1] - x_values[n-2]));
	}
	(*y2order)[n - 1] = (temp[n - 1] - (*y2order)[n - 1] * temp[n - 2]) / ((*y2order)[n - 1] * (*y2order)[n - 2] + 1.0);
	for (int i = n - 2; i >= 0; --i) {
		(*y2order)[i] = (*y2order)[i] * (*y2order)[i + 1] + temp[i];
	}
	return y2order;
}

void spline_solve_y(const vector<double> &x_values, const vector<double> &y_values, const vector<double> &y2order, const double x, double &y){
	int n = x_values.size();
	int x_boundary1 = 0,x_boundary2=n-1,pos;
	while (x_boundary2 - x_boundary1 > 1) {
		pos = (x_boundary1 + x_boundary2) / 2;
		if (x_values[pos] > x)x_boundary2 = pos;
		else x_boundary1 = pos;
	}
	double dis = x_values[x_boundary2] - x_values[x_boundary1];
	if (fabs(dis) < EPS) {
		cout << "��ֵ���󣺴�����ͬx" << endl;
		return;
	}
	double factor1 = (x_values[x_boundary2] - x) / dis;
	double factor2 = (x - x_values[x_boundary1]) / dis;
	y = factor1*y_values[x_boundary1] + factor2*y_values[x_boundary2] + ((factor1*factor1*factor1 - factor1)*y2order[x_boundary1] + (factor2*factor2*factor2 - factor2)*y2order[x_boundary2] * dis*dis) / 6.0;
}

