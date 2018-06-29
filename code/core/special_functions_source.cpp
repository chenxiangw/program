#include<iostream>
#include <cmath>
#include<vector>
#include<algorithm>
#include"Evaluation_function.h"
#define PI 3.141592653589793
#define EPS 1.0e-10
using namespace std;

bool Chebyshev_solve_factors(const double a, const double b, vector<double> &c, Evaluation_function &func) {
	/*
	计算切比雪夫多项式的系数
	参数：
	a：区间下界
	b：区间上界
	c：计算的切比雪夫多项式系数
	func：拟合函数
	*/
	double bpa, bma;
	int n = c.size();
	if (n == 0) {
		cout << "Chebyshev_solve_factors参数错误" << endl;
		return false;
	}
	vector<double>f = vector<double>(n);
	bma = 0.5*(b - a);
	bpa = 0.5*(b + a);
	for (int k = 0; k<n; ++k) {
		double y = cos(PI*(k + 0.5) / n);
		f[k] = func(y*bma + bpa);
	}
	double factor = 2.0 / n;
	for (int j = 0; j<n; ++j) {
		double sum = 0.0;
		for (int k = 0; k<n; ++k)sum += f[k] * cos(PI*j*(k + 0.5) / n);
		c[j] = factor * sum;
	}
	return true;
}

bool Chebyshev_value(const double a, const double b, const vector<double> &c, const double x, double &y) {
	/*
	计算切比雪夫多项式的值
	参数：
	a：区间下界
	b：区间上界
	c：切比雪夫多项式系数
	x：函数求值的位置
	*/
	int n = c.size();
	if ((x - a)*(x - b) > 0.0) {
		cout << "切比雪夫多项式求值错误：求值位置不在拟合区间内" << endl;
		return false;
	}
	double y2 = 2.0*(y = (2.0*x - a - b) / (b - a));
	double d = 0.0, pre_d = 0.0, temp;
	for (int j = n - 1; j >= 1; --j) {
		temp = d;
		d = y2 * d - pre_d + c[j];
		pre_d = temp;
	}
	y = y * d - pre_d + 0.5*c[0];
	return true;
}

double ln_gamma_Lanczos(const double x){
	/*
	Lanczos方法求gamma函数ln值
	Γ(z)=∫0->∞ t^(z-1)e^-tdt
	输入：
	x>0
	返回：
	ln(Γ(x))
	*/
	if (x <= 0) {
		cout << "ln_gamma_Lanczos参数错误" << endl;
		return false;
	}
	const int N = 6;//N=6时误差项<2*10^-10
	static double coefficient[N] = { 76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5 };

	double y = x;
	double tmp = x + 5.5;
	tmp -= (x + 0.5)*log(tmp);
	double value = 1.000000000190015;
	for (int j = 0; j < 6; ++j) value += coefficient[j] / ++y;
	return -tmp + log(sqrt(2*PI)*value / x);
}

double factorial(const int n){
	/*
	计算阶乘
	输入：
		n
	返回：
		n!
	*/
	static int ntop = 4;
	static double a[33] = { 1.0,1.0,2.0,6.0,24.0 };

	if (n < 0) cout << "阶乘参数错误" << endl;
	if (n > 32) return exp(ln_gamma_Lanczos(n + 1.0));//n>32时用gamma函数计算阶乘
	//n<=32时由表中计算
	int j;
	while (ntop<n) {
		j = ntop++;
		a[ntop] = a[j] * ntop;
	}
	return a[n];
}

double ln_factorial(const int n){
	/*
	计算阶乘的ln值
	输入：
	n
	返回：
	ln(n!)
	*/
	static double a[101];
	if (n < 0) cout << "阶乘参数错误" << endl;
	if (n <= 1) return 0.0;
	if (n <= 100) return a[n] ? a[n] : (a[n] = ln_gamma_Lanczos(n + 1.0));
	else return ln_gamma_Lanczos(n + 1.0);
}

double binomial(const int n, const int k){
	/*
	计算二项式系数Cnk
	输入：
	n,k
	返回：
	Cnk
	*/
	return floor(0.5 + exp(ln_factorial(n) - ln_factorial(k) - ln_factorial(n - k)));
}

double get_gamma_P(const double a, const double x){
	/*
	级数展开求求完全gamma函数P(a,x)
	输入：
	a,x
	输出：
	P：P(a,x)
	*/
	const int ITER_MAX = 100;
	double ln_gamma_a = ln_gamma_Lanczos(a);
	if (x < 0.0) {
		cout << "get_gamma_P参数错误" << endl;
		return 0.0;
	}
	double ap = a;
	double term = 1.0 / a;
	double sum = term;
	for (int n = 1; n <= ITER_MAX; ++n) {
		++ap;
		term *= x / ap;
		sum += term;
		if (fabs(term) < fabs(sum)*EPS) {
			return sum * exp(-x + a * log(x) - ln_gamma_a);
		}
	}
	cout << "get_gamma_P迭代未收敛" << endl;
	return sum * exp(-x + a * log(x) - ln_gamma_a);
}

double get_gamma_Q(const double a, const double x){
	/*
	连分式展开求Q(a,x)
	输入：
	a,x
	输出：
	Q：Q(a,x)
	*/
	const int ITER_MAX = 100;
	double ln_gamma_a = ln_gamma_Lanczos(a);
	if (x < 0.0) {
		cout << "get_gamma_Q参数错误" << endl;
		return 0.0;
	}
	double an, b, c, d, del, h;
	b = x + 1.0 - a;
	c = 1.0 / DBL_MIN;
	d = 1.0 / b;
	h = d;
	for (int i = 1; i <= ITER_MAX; ++i) {
		an = -i * (i - a);
		b += 2.0;
		d = an * d + b;
		if (fabs(d) < DBL_MIN) d = DBL_MIN;
		c = b + an / c;
		if (fabs(c) < DBL_MIN) c = DBL_MIN;
		d = 1.0 / d;
		del = d * c;
		h *= del;
		if (fabs(del - 1.0) < EPS) {
			return exp(-x + a * log(x) - ln_gamma_a)*h;
		}
	}
	cout << "get_gamma_Q迭代未收敛" << endl;
	return exp(-x + a * log(x) - ln_gamma_a)*h;
}

double gamma_P(const double a, const double x){
	/*
	求不完全gamma函数P(a,x)
	P(a,x)=γ(a,x)/Γ(a)=1/Γ(a)*∫0->x e^-t*t^(a-1)dt(a>0)
	输入：
	a,x
	返回：
	函数值
	*/

	if (x < 0.0 || a <= 0.0)  cout << "gamma_P参数错误" << endl;
	if (x < (a + 1.0)) {//用级数展开计算
		return get_gamma_P(a, x);
	}
	else {//用连分展开计算
		return 1.0 - get_gamma_Q(a, x);
	}
}

double gamma_Q(const double a, const double x){
	/*
	求不完全gamma函数P(a,x)的互补Q(a,x)
	Q(a,x)=1-P(a,x)
	输入：
	a,x
	返回：
	函数值
	*/

	if (x < 0.0 || a <= 0.0)  cout << "gamma_Q参数错误" << endl;
	if (x < (a + 1.0)) {//用级数展开计算
		return 1.0 - get_gamma_P(a, x);
	}
	else {//用连分展开计算
		return get_gamma_Q(a, x);
	}
}

double error_function(const double x){
	/*
	误差函数error_function(x)=P(1/2,x^2)(x>=0)
	输入：
	x
	返回：
	函数值
	*/
	return x < 0.0 ? -gamma_P(0.5, x*x) : gamma_P(0.5, x*x);
}

double exponential_integral(const int n, const double x){
	/*
	求指数积分函数En(x)
	En(x)=∫1->∞ (e^-xt)/t^n dt
	En(x)=x^(n-1)*Γ(1-n,x)
	输入：
	n,x
	返回：
	函数值
	*/
	const double EULER = 0.5772156649;
	const int ITER_MAX = 100;
	double result = 0.0;
	if (n < 0 || x < 0.0 || (x == 0.0 && (n == 0 || n == 1))){
		cout << "exponential_integral参数错误" << endl;
		return result;
	}
	if (n == 0) result = exp(-x) / x;
	else {
		if (fabs(x) < EPS) result = 1.0 / (n - 1);
		else {
			if (x > 1.0) {//用连分式计算
				double a, b, c, d, del, h;
				b = x + n;
				c = 1.0 / DBL_MIN;
				d = 1.0 / b;
				h = d;
				for (int i = 1; i <= ITER_MAX; ++i) {
					a = -i * (n - 1 + i);
					b += 2.0;
					d = 1.0 / (a*d + b);
					c = b + a / c;
					del = c * d;
					h *= del;
					if (fabs(del - 1.0) < EPS) {
						result = h * exp(-x);
						return result;
					}
				}
				cout << "exponential_integral迭代未收敛" << endl;
				return result;
			}
			else {//用幂级数计算
				result = (n - 1 != 0 ? 1.0 / (n - 1) : -log(x) - EULER);
				double factor = 1.0;
				double del, psi;
				for (int i = 1; i <= ITER_MAX; ++i) {
					factor *= -x / i;
					if (i != n - 1) del = -factor / (i - n + 1);
					else {
						psi = -EULER;
						for (int j = 1; j <= n - 1; ++j) psi += 1.0 / j;
						del = factor * (-log(x) + psi);
					}
					result += del;
					if (fabs(del) < fabs(result)*EPS) return result;
				}
				cout << "exponential_integral迭代未收敛" << endl;
				return result;
			}
		}
	}
	return result;
}

double exponential_integral_Ei(const double x){
	/*
	求Ei(x)
	Ei(x)=∫-∞->x (e^t)/t dt
	输入：
	x
	返回：
	函数值
	*/
	const int ITER_MAX = 100;
	const double EULER = 0.5772156649;
	if (x <= 0.0) {
		cout << "exponential_integral_Ei参数错误" << endl;
		return 0;
	}
	if (x < DBL_MIN) return log(x) + EULER;
	if (x <= -log(EPS)) {//幂级数求值
		double sum = 0.0;
		double factor = 1.0;
		double term;
		for (int k = 1; k <= ITER_MAX; ++k) {
			factor *= x / k;
			term = factor / k;
			sum += term;
			if (term < EPS*sum) {
				return sum + log(x) + EULER;
			}
		}
		cout << "exponential_integral_Ei迭代未收敛" << endl;
		return sum + log(x) + EULER;
	}
	else {//渐近级数求值
		double sum = 0.0;
		double term = 1.0;
		double prev;
		for (int k = 1; k <= ITER_MAX; ++k) {
			prev = term;
			term *= k / x;
			if (term < EPS) break;
			if (term < prev) sum += term;
			else {
				sum -= prev;
				break;
			}
		}
		return exp(x)*(1.0 + sum) / x;
	}
}

double beta(const double z, const double w) {
	/*
	求beta函数
	β(z,w)=β(w,z)=∫0->1 t^(z-1)(1-t)^(w-1)dt
	输入：
	z,w
	返回：
	函数值
	*/
	return exp(ln_gamma_Lanczos(z) + ln_gamma_Lanczos(w) - ln_gamma_Lanczos(z + w));
}

double incomplete_beta_continued_fractions(const double a, const double b, const double x){
	/*
	求不完全beta函数Ix(a,b)的连分式的值
	输入：
	a,b,x
	返回：
	连分式值
	*/
	const int ITER_MAX = 100;

	double an, c, d, del, h;
	c = 1.0;
	d = 1.0 - (a + b) * x / (a + 1.0);
	if (fabs(d) < DBL_MIN) d = DBL_MIN;
	d = 1.0 / d;
	h = d;
	for (int m = 1; m <= ITER_MAX; ++m) {
		double m2 = 2 * m;
		an = m * (b - m)*x / ((a - 1.0 + m2)*(a + m2));
		d = 1.0 + an * d;
		if (fabs(d) < DBL_MIN) d = DBL_MIN;
		c = 1.0 + an / c;
		if (fabs(c) < DBL_MIN) c = DBL_MIN;
		d = 1.0 / d;
		h *= d * c;
		an = -(a + m)*(a + b + m)*x / ((a + m2)*(a + 1.0 + m2));
		d = 1.0 + an * d;
		if (fabs(d) < DBL_MIN) d = DBL_MIN;
		c = 1.0 + an / c;
		if (fabs(c) < DBL_MIN) c = DBL_MIN;
		d = 1.0 / d;
		del = d * c;
		h *= del;
		if (fabs(del - 1.0) < EPS) return h;
	}
	cout << "exponential_integral_Ei迭代未收敛" << endl;
	return h;
}

double incomplete_beta(const double a, const double b, const double x){
	/*
	求不完全beta函数Ix(a,b)
	Ix(a,b)=∫0->x t^(a-1)*(1-t)^(b-1)dt/β(a,b)
	输入：
	a,b,x
	返回：
	函数值
	*/

	if (x < 0.0 || x > 1.0) {
		cout << "incomplete_beta参数错误" << endl;
		return 0.0;
	}
	double factor;
	if (x == 0.0 || x == 1.0) factor = 0.0;
	else factor = exp(ln_gamma_Lanczos(a + b) - ln_gamma_Lanczos(a) - ln_gamma_Lanczos(b) + a * log(x) + b * log(1.0 - x));
	if (x < (a + 1.0) / (a + b + 2.0))
		return factor * incomplete_beta_continued_fractions(a, b, x) / a;
	else
		return 1.0 - factor * incomplete_beta_continued_fractions(b, a, 1.0 - x) / b;
}

double Bessel_function_J0(const double x){
	/*
	求0阶一类贝塞尔函数J0(x)
	输入：
	x
	返回：
	函数值
	*/
	double y, result, factor1, factor2;
	if (fabs(x) < 8.0) {//有理函数逼近
		y = x * x;
		factor1 = 57568490574.0 + y * (-13362590354.0 + y * (651619640.7
			+ y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))));
		factor2 = 57568490411.0 + y * (1029532985.0 + y * (9494680.718
			+ y * (59272.64853 + y * (267.8532712 + y * 1.0))));
		result = factor1 / factor2;
	}
	else {//特定逼近函数
		y = 64.0 / x / x;
		factor1 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4
			+ y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
		factor2 = -0.1562499995e-1 + y * (0.1430488765e-3
			+ y * (-0.6911147651e-5 + y * (0.7621095161e-6
				- y * 0.934935152e-7)));
		result = sqrt(2 / PI / fabs(x))*(cos(fabs(x) - PI / 4)*factor1 - (8.0 / fabs(x)) * sin(fabs(x) - PI / 4)*factor2);
	}
	return result;
}

double Bessel_function_J1(const double x){
	/*
	求1阶一类贝塞尔函数J1(x)
	输入：
	x
	返回：
	函数值
	*/

	double y, result, factor1, factor2;
	if (fabs(x) < 8.0) {//有理函数逼近
		y = x * x;
		factor1 = x * (72362614232.0 + y * (-7895059235.0 + y * (242396853.1
			+ y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606))))));
		factor2 = 144725228442.0 + y * (2300535178.0 + y * (18583304.74
			+ y * (99447.43394 + y * (376.9991397 + y * 1.0))));
		result = factor1 / factor2;
	}
	else {//特定逼近函数
		y = 64.0 / x / x;
		factor1 = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4
			+ y * (0.2457520174e-5 + y * (-0.240337019e-6))));
		factor2 = 0.04687499995 + y * (-0.2002690873e-3
			+ y * (0.8449199096e-5 + y * (-0.88228987e-6
				+ y * 0.105787412e-6)));
		result = sqrt(2 / PI / fabs(x))*(cos(fabs(x) - PI * 3 / 4)*factor1 - (8.0 / fabs(x)) * sin(fabs(x) - PI * 3 / 4)*factor2);
		if (x < 0.0) result = -result;
	}
	return result;
}

double Bessel_function_J(const int n, const double x){
	/*
	求一类贝塞尔函数Jn(x)
	输入：
	n>=2,x
	返回：
	函数值
	*/
	if (n < 2 ) {
		cout << "Bessel_function_J参数错误" << endl;
		return 0.0;
	}
	double result;
	double abs_x = fabs(x);
	if (abs_x < EPS)return 0.0;
	else if (abs_x >(double) n) {//从0项向前递推
		double pre_J = Bessel_function_J0(abs_x);
		result = Bessel_function_J1(abs_x);
		double post_J;
		for (int j = 1; j<n; ++j) {
			post_J = j * 2.0 / abs_x *result - pre_J;
			pre_J = result;
			result = post_J;
		}
	}
	else {//从m项向后递推
		const double ACC = 40.0;//增加ACC提高精度
		const double BIGNO = 1.0e10;
		const double BIGNI = 1.0e-10;
		int m = 2 * ((n + (int)sqrt(ACC*n)) / 2);//m从n值大一个增量(常数*n)^1/2开始，其中常数的平方根粗略等于精度有效位数
		bool jsum = false;
		double sum = 0.0;
		double pre_J = result = 0.0;
		double post_J;
		double bJ = 1.0;
		for (int j = m; j>0; --j) {
			post_J = j * 2.0 / abs_x * bJ - pre_J;
			pre_J = bJ;
			bJ = post_J;
			if (fabs(bJ) > BIGNO) {//归一化防止溢出
				bJ *= BIGNI;
				pre_J *= BIGNI;
				result *= BIGNI;
				sum *= BIGNI;
			}
			if (jsum) sum += bJ;
			jsum = !jsum;
			if (j == n) result = pre_J;
		}
		sum = 2.0*sum - bJ;
		result /= sum;
	}
	return x < 0.0 && (n & 1) ? -result : result;
}

double Bessel_function_Y0(const double x){
	/*
	求0阶二类贝塞尔函数Y0(x)
	输入：
	x>0
	返回：
	函数值
	*/
	if (x <= 0) {
		cout << "Bessel_function_Y0参数错误" << endl;
		return 0.0;
	}

	double y, result, factor1, factor2;
	if (x < 8.0) {//有理函数逼近
		y = x * x;
		factor1 = -2957821389.0 + y * (7062834065.0 + y * (-512359803.6
			+ y * (10879881.29 + y * (-86327.92757 + y * 228.4622733))));
		factor2 = 40076544269.0 + y * (745249964.8 + y * (7189466.438
			+ y * (47447.26470 + y * (226.1030244 + y * 1.0))));
		result = (factor1 / factor2) + 0.636619772*Bessel_function_J0(x)*log(x);
	}
	else {//特定逼近函数
		y = (8.0 / x) * (8.0 / x);
		factor1 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4
			+ y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
		factor2 = -0.1562499995e-1 + y * (0.1430488765e-3
			+ y * (-0.6911147651e-5 + y * (0.7621095161e-6
				+ y * (-0.934945152e-7))));
		result = sqrt(2 / PI / x)*(sin(x  - PI / 4)*factor1 + (8.0 / x) * cos(x - PI / 4)*factor2);
	}
	return result;
}

double Bessel_function_Y1(const double x){
	/*
	求1阶二类贝塞尔函数Y1(x)
	输入：
	x>0
	返回：
	函数值
	*/
	if (x <= 0) {
		cout << "Bessel_function_Y1参数错误" << endl;
		return 0.0;
	}

	double y, result, factor1, factor2;
	if (x < 8.0) {//有理函数逼近
		y = x * x;
		factor1 = x * (-0.4900604943e13 + y * (0.1275274390e13
			+ y * (-0.5153438139e11 + y * (0.7349264551e9
				+ y * (-0.4237922726e7 + y * 0.8511937935e4)))));
		factor2 = 0.2499580570e14 + y * (0.4244419664e12
			+ y * (0.3733650367e10 + y * (0.2245904002e8
				+ y * (0.1020426050e6 + y * (0.3549632885e3 + y)))));
		result = (factor1 / factor2) + 0.636619772*(Bessel_function_J1(x)*log(x) - 1.0 / x);
	}
	else {//特定逼近函数
		y = (8.0 / x)  * (8.0 / x);
		factor1 = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4
			+ y * (0.2457520174e-5 + y * (-0.240337019e-6))));
		factor2 = 0.04687499995 + y * (-0.2002690873e-3
			+ y * (0.8449199096e-5 + y * (-0.88228987e-6
				+ y * 0.105787412e-6)));
		result = sqrt(2 / PI / x)*(sin(x - PI * 3 / 4)*factor1 + (8.0 / x)  * cos(x - PI * 3 / 4)*factor2);
	}
	return result;
}

double Bessel_function_Y(const int n, const double x){
	/*
	求二类贝塞尔函数Yn(x)
	输入：
	n>=2,x>0
	返回：
	函数值
	*/
	if (n < 2||x<=0) {
		cout << "Bessel_function_Y参数错误" << endl;
		return 0.0;
	}
	double result = Bessel_function_Y1(x);
	double pre_Y = Bessel_function_Y0(x);
	double post_Y;
	for (int j = 1; j<n; ++j) {
		post_Y = j * 2.0 / x * result - pre_Y;
		pre_Y = result;
		result = post_Y;
	}
	return result;
}

void Bessel_Chebyshev(const double x, double &gam1, double &gam2, double &gampl, double &gammi){
	/*
	|x|<=1/2时用切比雪夫展开式计算Γ1(x),Γ2(x)
	输入：
	x：求值位置
	输出：
	gam1：Γ1(x)
	gam2：Γ2(x)
	gampl：1/(Γ(1+x))
	gammi：1/(Γ(1-x))
	*/
	vector<double> c1 = {
		-1.142022680371172e0,6.516511267076e-3,
		3.08709017308e-4,-3.470626964e-6,6.943764e-9,
		3.6780e-11,-1.36e-13 };
	vector<double> c2 = {
		1.843740587300906e0,-0.076852840844786e0,
		1.271927136655e-3,-4.971736704e-6,-3.3126120e-8,
		2.42310e-10,-1.70e-13,-1.0e-15 };

	double new_x = 8.0*x*x - 1.0;
	Chebyshev_value(-1.0, 1.0, c1, new_x, gam1);
	Chebyshev_value(-1.0, 1.0, c2, new_x, gam2);
	gampl = gam2 - x * gam1;
	gammi = gam2 + x * gam1;
}

bool Bessel_function_JY(const double x, const double nu, double &output_J_nu, double &output_Y_nu, double &output_J_nu_diff, double &output_Y_nu_diff) {
	/*
	求分数阶一、二类贝塞尔函数Jν(x)和Yν(x)及其导数
	输入：
	x
	nu：ν
	输出：
	output_J_nu：Jν(x)
	output_Y_nu：Yν(x)
	output_J_nu_diff：Jν(x)'
	output_Y_nu_diff：Yν(x)'
	返回：
	计算无错误
	*/

	if (x <= 0.0 || nu < 0.0) {
		cout << "Bessel_function_JY参数错误" << endl;
		return false;
	}
	const double XMIN = 2.0;
	const int ITER_MAX = 10000;
	//x较小时采用Temme方法，只在|μ|<=1/2时有效，x较大时用Steed方法，保证μ<x即可
	int recursive_num = (x < XMIN ? (int)(nu + 0.5) : max(0, (int)(nu - x + 1.5)));//递推次数=ν-μ
	double mu = nu - recursive_num;//从ν递推到转折点μ<=x时收敛加快
	double xi = 1.0 / x;
	double w = 2.0*xi / PI;//Wronskian关系式
	//求第一个连分式fν=Jν'/Jν=...
	int isign = 1;
	double h = nu * xi;
	if (h < DBL_MIN) h = DBL_MIN;
	double b = 2.0*xi * nu;
	double d = 0.0;
	double c = h;
	double del;
	int i = 1;
	for (; i <= ITER_MAX; ++i) {
		b += 2.0*xi;
		d = b - d;
		if (fabs(d) < DBL_MIN) d = DBL_MIN;
		c = b - 1.0 / c;
		if (fabs(c) < DBL_MIN) c = DBL_MIN;
		d = 1.0 / d;
		del = c * d;
		h = del * h;
		if (d < 0.0) isign = -isign;
		if (fabs(del - 1.0) < EPS) break;
	}
	//h=fν
	if (i > ITER_MAX) {
		cout << "Bessel_function_JY迭代次数过多未收敛" << endl;
		return false;
	}
	//向后递推recursive_num项计算复连分式p+iq=(Jν'+iYν')/(Jν+iYν)=...
	double J_nu = isign * DBL_MIN;//向后递推时Jν先取任意值
	double J_nu_diff = h * J_nu;//Jν'=fν*Jν
	output_J_nu = J_nu;
	output_J_nu_diff = J_nu_diff;
	double factor = nu * xi;
	for (i = recursive_num; i >= 1; --i) {//Jν,Jν'递推到Jμ,Jμ'
		double pre_J_nu = factor * J_nu + J_nu_diff;
		factor -= xi;
		J_nu_diff = factor * pre_J_nu - J_nu;
		J_nu = pre_J_nu;
	}
	if (J_nu == 0.0) {
		cout << "Bessel_function_JY分母为0" << endl;
		return false;
	}
	double ratio = J_nu_diff / J_nu;//ratio=fμ

	double J_mu,Y_mu, Y_mu_diff,post_Y_mu, pre_Y_mu,p,q;
	if (x < XMIN) {//x较小时采用Temme方法
		double pimu = PI * mu;
		factor = (fabs(pimu) < EPS ? 1.0 : pimu / sin(pimu));
		d = -log(0.5*x);
		double e = mu * d;
		double gam1, gam2, gampl, gammi;
		Bessel_Chebyshev(mu, gam1, gam2, gampl, gammi);
		double ff = 2.0 / PI * factor*(gam1*cosh(e) + gam2 * (fabs(e) < EPS ? 1.0 : sinh(e) / e)*d);
		e = exp(e);
		p = e / (gampl*PI);
		q = 1.0 / (e*PI*gammi);
		factor = (fabs(0.5*pimu) < EPS ? 1.0 : sin(0.5*pimu) / (0.5*pimu));
		double r = PI * 0.5*pimu*factor*factor;
		c = 1.0;
		d = -0.25*x *x;
		double sum = ff + r * q;
		double sum1 = p;
		for (i = 1; i <= ITER_MAX; ++i) {
			ff = (i*ff + p + q) / (i*i - mu * mu);
			c *= (d / i);
			p /= (i - mu);
			q /= (i + mu);
			del = c * (ff + r * q);
			sum += del;
			sum1 += c * p - i * del;
			if (fabs(del) < (1.0 + fabs(sum))*EPS) break;
		}
		if (i > ITER_MAX) {
			cout << "Bessel_function_JY迭代次数过多未收敛" << endl;
			return false;
		}
		Y_mu = -sum;
		post_Y_mu = -sum1 * 2.0*xi;
		Y_mu_diff = mu * xi*Y_mu - post_Y_mu;
		J_mu = w / (Y_mu_diff - ratio * Y_mu);
	}
	else {//x较大时用Steed方法
		double a = 0.25 - mu * mu;
		p = -0.5*xi;
		q = 1.0;
		double br = 2.0*x;
		double bi = 2.0;
		factor = a * xi / (p*p + q * q);
		double cr = br + q * factor;
		double ci = bi + p * factor;
		double den = br * br + bi * bi;
		double dr = br / den;
		double di = -bi / den;
		double dlr = cr * dr - ci * di;
		double dli = cr * di + ci * dr;
		double temp = p * dlr - q * dli;
		q = p * dli + q * dlr;
		p = temp;
		for (i = 2; i <= ITER_MAX; ++i) {
			a += 2 * (i - 1);
			bi += 2.0;
			dr = a * dr + br;
			di = a * di + bi;
			if (fabs(dr) + fabs(di) < DBL_MIN) dr = DBL_MIN;
			factor = a / (cr*cr + ci * ci);
			cr = br + cr * factor;
			ci = bi - ci * factor;
			if (fabs(cr) + fabs(ci) < DBL_MIN) cr = DBL_MIN;
			den = dr * dr + di * di;
			dr /= den;
			di /= -den;
			dlr = cr * dr - ci * di;
			dli = cr * di + ci * dr;
			temp = p * dlr - q * dli;
			q = p * dli + q * dlr;
			p = temp;
			if (fabs(dlr - 1.0) + fabs(dli) < EPS) break;
		}
		if (i > ITER_MAX) {
			cout << "Bessel_function_JY迭代次数过多未收敛" << endl;
			return false;
		}
		double gamma = (p - ratio) / q;
		J_mu = sqrt(w / ((p - ratio)*gamma + q));
		if (J_nu < 0) J_mu = -J_mu;
		Y_mu = J_mu * gamma;
		Y_mu_diff = Y_mu * (p + q / gamma);
		post_Y_mu = mu * xi*Y_mu - Y_mu_diff;
	}

	factor = J_mu / J_nu;
	output_J_nu = output_J_nu * factor;
	output_J_nu_diff = output_J_nu_diff * factor;
	pre_Y_mu = Y_mu;
	Y_mu = post_Y_mu;
	for (i = 1; i <= recursive_num; ++i) {
		post_Y_mu = (mu + i)*2.0*xi*Y_mu - pre_Y_mu;
		pre_Y_mu = Y_mu;
		Y_mu = post_Y_mu;
	}
	output_Y_nu = pre_Y_mu;
	output_Y_nu_diff = nu * xi*pre_Y_mu - Y_mu;
	return true;
}

double Bessel_function_I0(const double x){
	/*
	求0阶修正一类贝塞尔函数I0(x)
	输入：
	x
	返回：
	函数值
	*/
	double result, y;
	double abs_x = fabs(x);

	if (abs_x < 3.75) {
		y = x / 3.75;
		y *= y;
		result = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492
			+ y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))));
	}
	else {
		y = 3.75 / abs_x;
		result = (exp(abs_x) / sqrt(abs_x))*(0.39894228 + y * (0.1328592e-1
			+ y * (0.225319e-2 + y * (-0.157565e-2 + y * (0.916281e-2
				+ y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1
					+ y * 0.392377e-2))))))));
	}
	return result;
}

double Bessel_function_I1(const double x){
	/*
	求1阶修正一类贝塞尔函数I1(x)
	输入：
	x
	返回：
	函数值
	*/
	double result, y;
	double abs_x = fabs(x);

	if (abs_x < 3.75) {
		y = x / 3.75;
		y *= y;
		result = abs_x * (0.5 + y * (0.87890594 + y * (0.51498869 + y * (0.15084934
			+ y * (0.2658733e-1 + y * (0.301532e-2 + y * 0.32411e-3))))));
	}
	else {
		y = 3.75 / abs_x;
		result = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1
			- y * 0.420059e-2));
		result = 0.39894228 + y * (-0.3988024e-1 + y * (-0.362018e-2
			+ y * (0.163801e-2 + y * (-0.1031555e-1 + y * result))));
		result *= (exp(abs_x) / sqrt(abs_x));
	}
	return x < 0.0 ? -result : result;
}

double Bessel_function_I(const int n, const double x){
	/*
	求修正一类贝塞尔函数In(x)
	输入：
	n>=2,x
	返回：
	函数值
	*/
	if (n < 2) {
		cout << "Bessel_function_I参数错误" << endl;
		return 0.0;
	}
	double result;
	double abs_x = fabs(x);
	if (abs_x < EPS)return 0.0;
	const double ACC = 40.0;//增加ACC提高精度
	const double BIGNO = 1.0e10;
	const double BIGNI = 1.0e-10;
	int m = 2 * ((n + (int)sqrt(ACC*n)) / 2);//m从n值大一个增量(常数*n)^1/2开始，其中常数的平方根粗略等于精度有效位数
	double pre_I = result = 0.0;
	double post_I;
	double bI = 1.0;

	for (int j = m; j>0; --j) {
		post_I = pre_I + j * 2.0 / abs_x *bI;
		pre_I = bI;
		bI = post_I;
		if (fabs(bI) > BIGNO) {
			result *= BIGNI;
			bI *= BIGNI;
			pre_I *= BIGNI;
		}
		if (j == n) result = pre_I;
	}
	result *= Bessel_function_I0(x) / bI;
	return x < 0.0 && (n & 1) ? -result : result;

}

double Bessel_function_K0(const double x){
	/*
	求0阶修正二类贝塞尔函数K0(x)
	输入：
	x>0
	返回：
	函数值
	*/

	double y, result;

	if (x <= 2.0) {
		y = x * x / 4.0;
		result = (-log(x / 2.0)*Bessel_function_I0(x)) + (-0.57721566 + y * (0.42278420
			+ y * (0.23069756 + y * (0.3488590e-1 + y * (0.262698e-2
				+ y * (0.10750e-3 + y * 0.74e-5))))));
	}
	else {
		y = 2.0 / x;
		result = (exp(-x) / sqrt(x))*(1.25331414 + y * (-0.7832358e-1
			+ y * (0.2189568e-1 + y * (-0.1062446e-1 + y * (0.587872e-2
				+ y * (-0.251540e-2 + y * 0.53208e-3))))));
	}
	return result;
}

double Bessel_function_K1(const double x){
	/*
	求1阶修正二类贝塞尔函数K1(x)
	输入：
	x>0
	返回：
	函数值
	*/

	double y, result;

	if (x <= 2.0) {
		y = x * x / 4.0;
		result = (log(x / 2.0)*Bessel_function_I1(x)) + (1.0 / x)*(1.0 + y * (0.15443144
			+ y * (-0.67278579 + y * (-0.18156897 + y * (-0.1919402e-1
				+ y * (-0.110404e-2 + y * (-0.4686e-4)))))));
	}
	else {
		y = 2.0 / x;
		result = (exp(-x) / sqrt(x))*(1.25331414 + y * (0.23498619
			+ y * (-0.3655620e-1 + y * (0.1504268e-1 + y * (-0.780353e-2
				+ y * (0.325614e-2 + y * (-0.68245e-3)))))));
	}
	return result;
}

double Bessel_function_K(const int n, const double x){
	/*
	求修正二类贝塞尔函数Kn(x)
	输入：
	n>=2,x>0
	返回：
	函数值
	*/
	if (n < 2 || x <= 0) {
		cout << "Bessel_function_K参数错误" << endl;
		return 0.0;
	}

	double result, pre_Y, post_Y;

	pre_Y = Bessel_function_K0(x);
	result = Bessel_function_K1(x);
	for (int j = 1; j<n; ++j) {
		post_Y = pre_Y + j * 2.0 / x * result;
		pre_Y = result;
		result = post_Y;
	}
	return result;
}



#undef PI
#undef EULER
#undef EPS