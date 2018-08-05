#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/legendre.hpp>
using namespace std;

double gamma(const double z) {
	/*
	Γ(z)=∫0->∞ t^(z-1)e^-tdt
	输入：
	z
	返回：
	Γ(z)
	*/
	return boost::math::tgamma(z);
}

double ln_gamma(const double z) {
	/*
	Γ(z)=∫0->∞ t^(z-1)e^-tdt
	输入：
	x
	返回：
	ln(|Γ(z)|)
	*/
	return boost::math::lgamma(z);
}

double gamma_P(const double a, const double z) {
	/*
	求不完全gamma函数P(a,x)
	P(a,z)=γ(a,z)/Γ(a)=1/Γ(a)*∫0->z e^-t*t^(a-1)dt(a>0)
	输入：
	a,z
	返回：
	P(a,z)
	*/
	if (a <= 0 || z < 0)  throw string("gamma_P参数错误");
	return boost::math::gamma_p(a, z);
}

double gamma_Q(const double a, const double z) {
	/*
	求不完全gamma函数P(a,x)
	P(a,z)=γ(a,z)/Γ(a)=1/Γ(a)*∫0->z e^-t*t^(a-1)dt(a>0)
	输入：
	a,z
	返回：
	P(a,z)
	*/
	if (a <= 0 || z < 0)  throw string("gamma_Q参数错误");
	return boost::math::gamma_q(a, z);
}

double factorial(const int n) {
	/*
	计算阶乘
	输入：
	n
	返回：
	n!
	*/
	return boost::math::factorial<double>(n);
}

double binomial(const int n, const int k) {
	/*
	计算二项式系数Cnk
	输入：
	n,k
	返回：
	Cnk
	*/
	return boost::math::binomial_coefficient<double>(n, k);
}

double beta(const double a, const double b) {
	/*
	求beta函数
	β(a,b)=Γ(a)Γ(b)/Γ(a+b)=∫0->1 t^(a-1)(1-t)^(b-1)dt
	输入：
	a,b
	返回：
	β(a,b)
	*/
	return boost::math::beta(a, b);
}

double incomplete_beta(const double a, const double b, const double x) {
	/*
	求不完全beta函数Ix(a,b)
	Ix(a,b)=∫0->x t^(a-1)*(1-t)^(b-1)dt/β(a,b)
	输入：
	a,b,x
	返回：
	Ix(a,b)
	*/
	return boost::math::beta(a, b, x);
}

double error_function(const double x) {
	/*
	erf(x)=P(1/2,x^2)(x>=0)
	输入：
	x
	返回：
	erf(x)
	*/
	return boost::math::erf(x);
}

double Legendre_p(const int l,const double x) {
	/*
	勒让德多项式
	输入：
	l,x
	返回：
	Pl(x)
	*/
	if (x < -1 || x > 1)  throw string("legendre_p参数错误");
	return boost::math::legendre_p(l,x);
}

double Legendre_q(const int l, const double x) {
	/*
	第二类勒让德多项式
	输入：
	l,x
	返回：
	Ql(x)
	*/
	if (x < -1 || x > 1)  throw string("legendre_p参数错误");
	return boost::math::legendre_q(l,x);
}

