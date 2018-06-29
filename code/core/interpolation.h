#pragma once
#include<vector>
using std::vector;

double polynomial_interpolation(const vector<double> &x_values, const vector<double> &y_values, const double x, double &y, const int start = 0, const int end = -1);

double rational_interpolation(const vector<double> &x_values, const vector<double> &y_values, const double x, double &y);

void spline_solve_y(const vector<double> &x_values, const vector<double> &y_values, const vector<double> &y2order, const double x, double &y);