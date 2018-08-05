#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/legendre.hpp>
using namespace std;

double gamma(const double z) {
	/*
	��(z)=��0->�� t^(z-1)e^-tdt
	���룺
	z
	���أ�
	��(z)
	*/
	return boost::math::tgamma(z);
}

double ln_gamma(const double z) {
	/*
	��(z)=��0->�� t^(z-1)e^-tdt
	���룺
	x
	���أ�
	ln(|��(z)|)
	*/
	return boost::math::lgamma(z);
}

double gamma_P(const double a, const double z) {
	/*
	����ȫgamma����P(a,x)
	P(a,z)=��(a,z)/��(a)=1/��(a)*��0->z e^-t*t^(a-1)dt(a>0)
	���룺
	a,z
	���أ�
	P(a,z)
	*/
	if (a <= 0 || z < 0)  throw string("gamma_P��������");
	return boost::math::gamma_p(a, z);
}

double gamma_Q(const double a, const double z) {
	/*
	����ȫgamma����P(a,x)
	P(a,z)=��(a,z)/��(a)=1/��(a)*��0->z e^-t*t^(a-1)dt(a>0)
	���룺
	a,z
	���أ�
	P(a,z)
	*/
	if (a <= 0 || z < 0)  throw string("gamma_Q��������");
	return boost::math::gamma_q(a, z);
}

double factorial(const int n) {
	/*
	����׳�
	���룺
	n
	���أ�
	n!
	*/
	return boost::math::factorial<double>(n);
}

double binomial(const int n, const int k) {
	/*
	�������ʽϵ��Cnk
	���룺
	n,k
	���أ�
	Cnk
	*/
	return boost::math::binomial_coefficient<double>(n, k);
}

double beta(const double a, const double b) {
	/*
	��beta����
	��(a,b)=��(a)��(b)/��(a+b)=��0->1 t^(a-1)(1-t)^(b-1)dt
	���룺
	a,b
	���أ�
	��(a,b)
	*/
	return boost::math::beta(a, b);
}

double incomplete_beta(const double a, const double b, const double x) {
	/*
	����ȫbeta����Ix(a,b)
	Ix(a,b)=��0->x t^(a-1)*(1-t)^(b-1)dt/��(a,b)
	���룺
	a,b,x
	���أ�
	Ix(a,b)
	*/
	return boost::math::beta(a, b, x);
}

double error_function(const double x) {
	/*
	erf(x)=P(1/2,x^2)(x>=0)
	���룺
	x
	���أ�
	erf(x)
	*/
	return boost::math::erf(x);
}

double Legendre_p(const int l,const double x) {
	/*
	���õ¶���ʽ
	���룺
	l,x
	���أ�
	Pl(x)
	*/
	if (x < -1 || x > 1)  throw string("legendre_p��������");
	return boost::math::legendre_p(l,x);
}

double Legendre_q(const int l, const double x) {
	/*
	�ڶ������õ¶���ʽ
	���룺
	l,x
	���أ�
	Ql(x)
	*/
	if (x < -1 || x > 1)  throw string("legendre_p��������");
	return boost::math::legendre_q(l,x);
}

