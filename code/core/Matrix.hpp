#pragma once
#include <vector>
#include<iostream>
#include<memory>
#include<algorithm>
#define EPS pow(2.0, -52.0)

template <typename value_type>
class Matrix
{
private:
	std::vector<std::vector<value_type>> *data = nullptr;

	inline std::shared_ptr<Matrix<double>> Matrix<value_type>::inverse_LU_method(double eps = pow(2.0, -52.0)) const;
	inline double Matrix<value_type>::determinant_LU_method() const;
public:
	//构造函数
	Matrix();//默认
	Matrix(int m,int n);//m*n零矩阵
	Matrix(std::vector<std::vector<value_type>> &input);//vector复制构造
	Matrix(const Matrix<value_type> &input);//Matrix复制构造
	//赋值
	Matrix<value_type> & operator=(std::vector<std::vector<value_type>> &input);//=vector
	Matrix<value_type> & operator=(const Matrix<value_type> &input);//=Matrix
	//[]重载
	inline std::vector<value_type> & operator[](const int i);
	inline const std::vector<value_type> & operator[](const int i) const;//只读
	//运算重载
	 std::shared_ptr<Matrix<value_type>> operator+(Matrix<value_type> &b);
	 std::shared_ptr<Matrix<value_type>>  operator-(Matrix<value_type> &b);
	 std::shared_ptr<Matrix<value_type>>  operator*(Matrix<value_type> &b);
	 std::shared_ptr<Matrix<value_type>>  operator/(Matrix<value_type> &b);
	 std::shared_ptr<Matrix<value_type>> operator+(double b);
	 std::shared_ptr<Matrix<value_type>>  operator-(double b);
	 std::shared_ptr<Matrix<value_type>>  operator*(double b);
	 std::shared_ptr<Matrix<value_type>>  operator/(double b);
	 std::shared_ptr<Matrix<value_type>> reverse_sign();
	 Matrix<value_type> & operator+=(Matrix<value_type> &b);
	 Matrix<value_type> &  operator-=(Matrix<value_type> &b);
	 Matrix<value_type>  & operator*=(Matrix<value_type> &b);
	 Matrix<value_type> & operator+=(double b);
	 Matrix<value_type> &  operator-=(double b);
	 Matrix<value_type> &  operator*=(double b);
	//
	inline int msize() const;
	inline int nsize() const;
	//打印
	 void print() const;
	//析构
	~Matrix();

	//转置
	std::shared_ptr<Matrix<value_type>> transposition()const;
	//求逆
	std::shared_ptr<Matrix<double>> Matrix<value_type>::inverse()const;
	//求行列式
	double determinant()const;
};

template <typename value_type>
Matrix<value_type>::Matrix() {
	data = nullptr;
}

template <typename value_type>
Matrix<value_type>::Matrix(int m, int n) {
	data = new std::vector<std::vector<value_type>>();
	for (int i = 0; i < m; ++i) {
		std::vector<value_type> line(n,0);
		(*data).push_back(line);
	}
}

template <typename value_type>
Matrix<value_type>::Matrix(std::vector<std::vector<value_type>> &input) {
	int msize = input.size();
	int nsize = 0;
	if (msize != 0) nsize = input[0].size();
	for (int i = 1; i < msize; ++i) {
		if (input[i].size() != nsize) {
			std::cout << "输入格式错误" << endl;
			return;
		}
	}
	data = new std::vector<std::vector<value_type>>();
	for (int i = 0; i < msize; ++i) {
		std::vector<value_type> line(nsize);
		for (int j = 0; j < nsize; ++j) {
			line[j] = input[i][j];
		}
		(*data).push_back(line);
	}
}

template <typename value_type>
Matrix<value_type>::Matrix(const Matrix<value_type> &input) {
	int msize = input.msize();
	int nsize = input.nsize();
	data = new std::vector<std::vector<value_type>>();
	for (int i = 0; i < msize; ++i) {
		std::vector<value_type>line(nsize);
		for (int j = 0; j < nsize; ++j) {
			line[j] = input[i][j];
		}
		(*data).push_back(line);
	}

}

template <typename value_type>
Matrix<value_type> & Matrix<value_type>::operator=(std::vector<std::vector<value_type>> &input) {
	int msize = input.size();
	int nsize = 0;
	if (msize != 0) nsize = input[0].size();
	for (int i = 1; i < msize; ++i) {
		if (input[i].size() != nsize) {
			std::cout << "输入格式错误" << std::endl;
			return *this;
		}
	}
	if (this->msize() != msize || this->nsize() != nsize) {
		if (data != nullptr) delete data;
		data = new std::vector<std::vector<value_type>>();
		for (int i = 0; i < msize; ++i) {
			std::vector<value_type> line(nsize);
			(*data).push_back(line);
		}
		for (int i = 0; i < msize; ++i) {
			for (int j = 0; j < nsize; ++j) {
				(*data)[i][j] = input[i][j];
			}
		}
	}
	return *this;
}

template <typename value_type>
Matrix<value_type> & Matrix<value_type>::operator=(const Matrix<value_type> &input) {
	int msize = input.msize();
	int nsize = input.nsize();
	if (this != &input) {
		if (this->msize() != msize|| this->nsize() != nsize) {
			if (data != nullptr) delete data;
			data = new std::vector<std::vector<value_type>>();
			for (int i = 0; i < msize; ++i) {
				std::vector<value_type> line(nsize);
				(*data).push_back(line);
			}
		}
		for (int i = 0; i < msize; ++i) {
			for (int j = 0; j < nsize; ++j) {
				(*data)[i][j] = input[i][j];
			}
		}
	}
	return *this;
}

template <typename value_type>
inline std::vector<value_type> & Matrix<value_type>::operator[](const int i) {
	return (*data)[i];
}

template <typename value_type>
inline const std::vector<value_type> & Matrix<value_type>::operator[](const int i) const {
	return  (*data)[i];
}

template <typename value_type>
 std::shared_ptr<Matrix<value_type>> Matrix<value_type>::operator+(Matrix<value_type> &b) {
	int m = this->msize();
	int n = this->nsize();
	std::shared_ptr<Matrix<value_type>> result(new Matrix<double>(m, n));
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			(*result)[i][j] = (*this)[i][j] + b[i][j];
		}
	}
	return result;
}
template <typename value_type>
 std::shared_ptr<Matrix<value_type>> Matrix<value_type>::operator-(Matrix<value_type> &b) {
	int m = this->msize();
	int n = this->nsize();
	std::shared_ptr<Matrix<value_type>> result(new Matrix<double>(m, n));
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			(*result)[i][j] = (*this)[i][j] - b[i][j];
		}
	}
	return result;
}
template <typename value_type>
 std::shared_ptr<Matrix<value_type>> Matrix<value_type>::operator*(Matrix<value_type> &b) {
	int m = this->msize();
	int n = b.nsize();
	int mid = this->nsize();
	std::shared_ptr<Matrix<value_type>> result(new Matrix<double>(m, n));
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			double sum = 0;
			for (int k = 0; k < mid; ++k) {
				sum += (*this)[i][k] * b[k][j];
			}
			(*result)[i][j] = sum;
		}
	}
	return result;
}
 template <typename value_type>
 std::shared_ptr<Matrix<value_type>> Matrix<value_type>::operator/(Matrix<value_type> &b) {
	 int m = this->msize();
	 int n = b.nsize();
	 int mid = this->nsize();
	 try {
		 std::shared_ptr<Matrix<value_type>> inverse = b.inverse();
		 std::shared_ptr<Matrix<value_type>> result(new Matrix<double>(m, n));
		 for (int i = 0; i < m; ++i) {
			 for (int j = 0; j < n; ++j) {
				 double sum = 0;
				 for (int k = 0; k < mid; ++k) {
					 sum += (*this)[i][k] * (*inverse)[k][j];
				 }
				 (*result)[i][j] = sum;
			 }
		 }
		 return result;
	 }
	 catch (std::string error) {
		 std::cout << error << std::endl;
		 return nullptr;
	 }
 }
 template <typename value_type>
 std::shared_ptr<Matrix<value_type>> Matrix<value_type>::operator+(double b) {
	 int m = this->msize();
	 int n = this->nsize();
	 std::shared_ptr<Matrix<value_type>> result(new Matrix<double>(m, n));
	 for (int i = 0; i < m; ++i) {
		 for (int j = 0; j < n; ++j) {
			 (*result)[i][j] = (*this)[i][j] + b;
		 }
	 }
	 return result;
 }
 template <typename value_type>
 std::shared_ptr<Matrix<value_type>>  Matrix<value_type>::operator-(double b) {
	 int m = this->msize();
	 int n = this->nsize();
	 std::shared_ptr<Matrix<value_type>> result(new Matrix<double>(m, n));
	 for (int i = 0; i < m; ++i) {
		 for (int j = 0; j < n; ++j) {
			 (*result)[i][j] = (*this)[i][j] - b;
		 }
	 }
	 return result;
 }
 template <typename value_type>
 std::shared_ptr<Matrix<value_type>> Matrix<value_type>:: operator*(double b) {
	 int m = this->msize();
	 int n = this->nsize();
	 std::shared_ptr<Matrix<value_type>> result(new Matrix<double>(m, n));
	 for (int i = 0; i < m; ++i) {
		 for (int j = 0; j < n; ++j) {
			 (*result)[i][j] = (*this)[i][j] * b;
		 }
	 }
	 return result;
 }
 template <typename value_type>
 std::shared_ptr<Matrix<value_type>> Matrix<value_type>:: operator/(double b) {
	 int m = this->msize();
	 int n = this->nsize();
	 std::shared_ptr<Matrix<value_type>> result(new Matrix<double>(m, n));
	 for (int i = 0; i < m; ++i) {
		 for (int j = 0; j < n; ++j) {
			 (*result)[i][j] = (*this)[i][j] / b;
		 }
	 }
	 return result;
 }
 template <typename value_type>
 std::shared_ptr<Matrix<value_type>>  Matrix<value_type>::reverse_sign() {
	 int m = this->msize();
	 int n = this->nsize();
	 std::shared_ptr<Matrix<value_type>> result(new Matrix<double>(m, n));
	 for (int i = 0; i < m; ++i) {
		 for (int j = 0; j < n; ++j) {
			 (*result)[i][j] = -(*this)[i][j];
		 }
	 }
	 return  result;
 }

 template <typename value_type>
 Matrix<value_type> & Matrix<value_type>::operator+=(Matrix<value_type> &b) {
	 int m = this->msize();
	 int n = this->nsize();
	 for (int i = 0; i < m; ++i) {
		 for (int j = 0; j < n; ++j) {
			 (*this)[i][j] += b[i][j];
		 }
	 }
	 return  (*this);
 }
 template <typename value_type>
 Matrix<value_type> & Matrix<value_type>::operator-=(Matrix<value_type> &b) {
	 int m = this->msize();
	 int n = this->nsize();
	 for (int i = 0; i < m; ++i) {
		 for (int j = 0; j < n; ++j) {
			 (*this)[i][j] -= b[i][j];
		 }
	 }
	 return  (*this);
 }
 template <typename value_type>
 Matrix<value_type> & Matrix<value_type>::operator*=(Matrix<value_type> &b) {
	 int m = this->msize();
	 int n = b.nsize();
	 int mid = this->nsize();
	 Matrix<value_type> result(m, n);
	 for (int i = 0; i < m; ++i) {
		 for (int j = 0; j < n; ++j) {
			 double sum = 0;
			 for (int k = 0; k < mid; ++k) {
				 sum += (*this)[i][k] * b[k][j];
			 }
			 result[i][j] = sum;
		 }
	 }
	 (*this) = result;
	 return (*this);
 }
 template <typename value_type>
 Matrix<value_type> & Matrix<value_type>::operator+=(double b) {
	 int m = this->msize();
	 int n = this->nsize();
	 for (int i = 0; i < m; ++i) {
		 for (int j = 0; j < n; ++j) {
			 (*this)[i][j] += b;
		 }
	 }
	 return (*this);
 }
 template <typename value_type>
 Matrix<value_type> &  Matrix<value_type>::operator-=(double b) {
	 int m = this->msize();
	 int n = this->nsize();
	 for (int i = 0; i < m; ++i) {
		 for (int j = 0; j < n; ++j) {
			 (*this)[i][j] -= b;
		 }
	 }
	 return (*this);
 }
 template <typename value_type>
 Matrix<value_type> &  Matrix<value_type>::operator*=(double b) {
	 int m = this->msize();
	 int n = this->nsize();
	 for (int i = 0; i < m; ++i) {
		 for (int j = 0; j < n; ++j) {
			 (*this)[i][j] *= b;
		 }
	 }
	 return (*this);
 }


template <typename value_type>
inline int Matrix<value_type>::msize() const {
	if (data == nullptr) return 0;
	return data->size();
}
template <typename value_type>
inline int Matrix<value_type>::nsize() const {
	if (data == nullptr) return 0;
	return (*data)[0].size();
}

template <typename value_type>
Matrix<value_type>::~Matrix() {
	if (data != nullptr) delete data;
}

template <typename value_type>
 void Matrix<value_type>::print() const {
	int msize = this->msize();
	int nsize = this->nsize();
	std::cout << "[";
	for (int i = 0; i < msize; ++i) {
		std::cout << "[ ";
		for (int j = 0; j < nsize; ++j) {
			std::cout << (*data)[i][j] << ' ';
		}
		std::cout << "]";
	}
	std::cout << "]" << std::endl;
}

template <typename value_type>
std::shared_ptr<Matrix<value_type>> Matrix<value_type>::transposition() const {
	std::shared_ptr<Matrix<value_type>> AT;
	if (this->msize() == 0||this->nsize()==0) {
		std::cout << "转置错误，矩阵为空" << endl;
		return AT;
	}

	int n = this->msize();
	int m = this->nsize();
	AT.reset(new Matrix<value_type>(m, n));
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			(*AT)[j][i] = (*this)[i][j];
		}
	}
	return AT;
}

bool LU_decomposition(Matrix<double> &A, std::vector<int> &index, double eps = EPS);
void LU_solve(const Matrix<double> &LU, const std::vector<int> &index, std::vector<double> &b);

template <typename value_type>
 std::shared_ptr<Matrix<double>> Matrix<value_type>::inverse_LU_method(double eps) const {

	int m = this->msize();
	std::vector<int> index(m);
	Matrix<double> LU(m,m);
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < m; ++j) {
			LU[i][j] = (double) (*data)[i][j];
		}
	}
	LU_decomposition(LU, index, eps);
	std::shared_ptr<Matrix<double>> inverse;
	double value = 1.0;
	for (int j = 0; j < m; ++j) {
		value *= LU[j][j];
	}
	if (value < eps) {
		std::cout << "求逆错误：矩阵为奇异矩阵" << std::endl;
		return inverse;
	}
	inverse.reset(new Matrix<double>(m, m));
	std::vector<double> temp_b(m);
	for (int j = 0; j < m; ++j) {
		for (int i = 0; i < m; ++i) temp_b[i] = 0.0;
		temp_b[j] = 1.0;
		LU_solve(LU, index, temp_b);
		for (int i = 0; i < m; ++i) (*inverse)[i][j] = temp_b[i];
	}
	return inverse;
}

template <typename value_type>
std::shared_ptr<Matrix<double>> Matrix<value_type>::inverse() const {
	std::shared_ptr<Matrix<double>> inverse;
	if (this->msize() == 0) {
		std::cout << "求逆错误：矩阵为空" << std::endl;
		return inverse;
	}
	if (this->msize() != this->nsize()) { 
		std::cout << "求逆错误：矩阵不是方阵" << std::endl;
		return inverse;
	}
	inverse= this->inverse_LU_method();
	return inverse;
}

template <typename value_type>
 double Matrix<value_type>::determinant_LU_method() const {
	/*
	求行列式
	参数：
	A：m*m系数矩阵
	返回：
	double：行列式值
	*/

	double value = 1.0;
	int m = this->msize();
	std::vector<int> index(m);
	Matrix<double> LU(m, m);
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < m; ++j) {
			LU[i][j] = (double)(*data)[i][j];
		}
	}
	bool even_exchange = LU_decomposition(LU, index);
	for (int j = 0; j < m; ++j) {
		value *= LU[j][j];
	}
	return even_exchange ? value : -value;
}

template <typename value_type>
double Matrix<value_type>::determinant() const {
	if (this->msize() == 0) {
		std::cout << "行列式错误：矩阵为空" << std::endl;
		return 0;
	}
	if (this->msize() != this->nsize()) {
		std::cout << "行列式错误：矩阵不是方阵" << std::endl;
		return 0;
	}
	return this->determinant_LU_method();
}

#undef EPS