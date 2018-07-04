#include<cmath>
#include<vector>
#include<string>
#include<algorithm>
#include<memory>
#include<deque>
#include"Matrix.hpp"
#include"Linear_solution.hpp"
#define EPS pow(10,-10)

using namespace std;

void Gauss_Jordan(Matrix<double> &A, Matrix<double> &b) {
	/*
		��ȫ��Ԫ����Gauss_Jordan��ȥ��
		������
			A��m*mϵ������
			b��m*n�Ҷ���
		���أ�
			Ax=b
			A������������
			b��x����
		ʱ�临�Ӷȣ�
			m^3+nm^2
	*/
	if (A.msize() == 0)throw std::string("ϵ������Ϊ��");
	if (A.msize() != A.nsize()) throw std::string("ϵ�������Ƿ���");

	size_t m = A.msize();
	size_t n = b.nsize();
	vector<int> index_col(m), index_row(m), pivoted(m);
	for (int i = 0; i < m; ++i) pivoted[i] = 0;
	for (int i = 0; i < m; ++i) {
		double max_num = 0.0;
		int pivot_row, pivot_col;
		//ѡȡA�����Ԫ��Ϊ��Ԫ����¼λ��
		for (int j = 0; j < m; ++j) {
			if (pivoted[j] != 1) {
				for (int k = 0; k < m; ++k) {
					if (pivoted[k] == 0) {
						if (fabs(A[j][k]) > max_num) {
							max_num = fabs(A[j][k]);
							pivot_row = j;
							pivot_col = k;
						}
					}
				}
			}
		}
		++(pivoted[pivot_col]);
		//����Ԫ���е��Խ���λ��
		if (pivot_row != pivot_col) {
			for (int k = 0; k < m; ++k)swap(A[pivot_row][k], A[pivot_col][k]);
			for (int k = 0; k < n; ++k)swap(b[pivot_row][k], b[pivot_col][k]);
		}
		index_row[i] = pivot_row;
		index_col[i] = pivot_col;
		//��Ԫ�����й�һ��
		if(A[pivot_col][pivot_col]==0.0) throw std::string("�������");
		double pivot_inverse = 1.0 / A[pivot_col][pivot_col];
		A[pivot_col][pivot_col] = 1.0;
		for (int k = 0; k < m; ++k)A[pivot_col][k] *= pivot_inverse;
		for (int k = 0; k < n; ++k)b[pivot_col][k] *= pivot_inverse;
		//����Ԫ����ȥ��Ԫ������
		for (int j = 0; j < m; ++j) {
			if (j != pivot_col) {
				double factor = A[j][pivot_col];
				A[j][pivot_col] = 0.0;
				for (int k = 0; k < m; ++k) A[j][k] -= A[pivot_col][k] * factor;
				for (int k = 0; k < n; ++k) b[j][k] -= b[pivot_col][k] * factor;
			}
		}
	}
	//��������д���ԭ
	for (int l = m - 1; l >= 0; --l) {
		if (index_row[l] != index_col[l]) {
			for (int k = 0; k < m; ++k) {
				swap(A[k][index_row[l]], A[k][index_col[l]]);
			}
		}
	}
}

bool LU_decomposition(Matrix<double> &A, vector<int> &index, double eps) {
	/*
		Crout�ֽⷨ
		������
		A��m*mϵ������
		index����¼���������
		���أ�
		LU=A
		A�����Ǳ�L����������Ǳ�U���
		bool���н��������Ƿ�Ϊż��
		ʱ�临�Ӷȣ�
		1/3m^3
	*/
	bool even_exchange = true;
	size_t m = A.msize();
	vector<double> factor_at_line(m);
	for (int i = 0; i < m; ++i) {
		double max_num = 0.0;
		for (int j = 0; j < m; ++j) {
			if (fabs(A[i][j]) > max_num) max_num = fabs(A[i][j]);
		}
		if (max_num == 0.0)  throw std::string("�������");
		factor_at_line[i] = 1.0 / max_num;
	}
	for (int j = 0; j < m; ++j) {
		//�����λ�õ�ֵ
		for (int i = 0; i < j; ++i) {
			double sum = A[i][j];
			for (int k = 0; k < i; ++k)sum -= A[i][k] * A[k][j];
			A[i][j] = sum;
		}
		//�����λ�õ�ֵ
		double max_num = 0.0;
		int max_pos;
		for (int i = j; i < m; ++i) {
			double sum = A[i][j];
			for (int k = 0; k < j; ++k) sum -= A[i][k] * A[k][j];
			A[i][j] = sum;
			double temp = factor_at_line[i] * fabs(sum);
			if (temp >= max_num) {
				max_num = temp;
				max_pos = i;
			}
		}
		if (j != max_pos) {
			for (int k = 0; k < m; ++k) {
				swap(A[max_pos][k], A[j][k]);
			}
			even_exchange = !even_exchange;
			swap(factor_at_line[max_pos], factor_at_line[j]);
		}
		index[j] = max_pos;
		if (A[j][j] < eps) {
			A[j][j] = 0;
			return true;
		}
		if (j != m - 1) {
			double factor = 1.0 / (A[j][j]);
			for (int i = j + 1; i < m; ++i) A[i][j] *= factor;
		}
	}
	return even_exchange;
}

void LU_solve(const Matrix<double> &LU,const vector<int> &index, vector<double> &b) {
	/*
		LU=b���
		������
		LU��LU���Ͼ���
		index����¼���������
		b���Ҷ�������
		���أ�
		b�������
		ʱ�临�Ӷȣ�
		2/3m^2
	*/
	double sum;
	int pos = 0;
	size_t m = LU.msize();
	for (int i = 0; i < m; ++i) {
		int ip = index[i];
		sum = b[ip];
		b[ip] = b[i];
		if (pos != 0)
			for (int j = pos - 1; j < i; j++) sum -= LU[i][j] * b[j];
		else if (sum != 0.0)
			pos = i + 1;
		b[i] = sum;
	}
	for (int i = m - 1; i >= 0; i--) {
		sum = b[i];
		for (int j = i + 1; j < m; j++) sum -= LU[i][j] * b[j];
		b[i] = sum / LU[i][i];
	}

}

void improve_LU_solve(const Matrix<double> &A, const Matrix<double> &LU, const vector<int> &index,const vector<double> &b, vector<double> &x){
	/*
		�Խ�x�Ż�
		������
		A��ԭϵ������
		LU��LU���Ͼ���
		index����¼���������
		b���Ҷ�������
		x��ԭ��
		���أ�
		b�������
	*/
	size_t m = A.msize();
	vector<double> differ(m);
	for (int i = 0; i < m; ++i) {
		long double temp_differ = -b[i];
		for (int j = 0; j < m; ++j) {
			temp_differ += (long double)A[i][j] * (long double)x[j];
		}
		differ[i] = temp_differ;
	}
	LU_solve(LU, index, differ);
	for (int i = 0; i < m; ++i) {
		x[i] -= differ[i];
	}
}

void solve_linear_LU_method(const Matrix<double> &A, Matrix<double> &b) {
	/*
		���ϵ�������������Է�����
		������
		A��m*mϵ������
		b��m*n�Ҷ���
		���أ�
		b����
	*/

	if (A.msize() == 0)throw "ϵ������Ϊ��";
	if (A.msize() != A.nsize()) throw "ϵ�������Ƿ���";
	if (b.nsize() == 0) throw "�Ҷ���Ϊ��";

	size_t m = A.msize();
	size_t n = b.nsize();
	vector<int> index(m);
	Matrix<double> LU(A);
	LU_decomposition(LU, index);//A�ֽ��LU����
	for (int i = 0; i < n; ++i) {
		vector<double> temp_b(m);
		vector<double> x(m);
		for (int j = 0; j < m; ++j) {
			x[j] = b[j][i];
			temp_b[j] = b[j][i];
		}
		LU_solve(LU, index, x);//���
		improve_LU_solve(A, LU, index, temp_b, x);//�Ż���
		for (int j = 0; j < m; ++j) {
			b[j][i] = x[j];
		}
	}

}

shared_ptr<Matrix<double>> matrix_mutiply(const Matrix<double> &A, const Matrix<double> &B) {
	if (A.nsize() != B.msize()) throw "�����޷����";
	size_t m = A.msize();
	size_t n = B.nsize();
	size_t mid = A.nsize();
	shared_ptr<Matrix<double>> AB(new Matrix<double>(m, n));
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			double sum = 0;
			for (int k = 0; k < mid; ++k) {
				sum += A[i][k] * B[k][j];
			}
			(*AB)[i][j] = sum;
		}
	}
	return AB;
}

void householder_transform(Matrix<double> &A, Matrix<double> &U, vector<double> &S, Matrix<double> &V, vector<double> &e) {
	/*
		householder��A�任Ϊ���ԽǾ���
		A=US'V������SΪS'�ĶԽ�Ԫ�أ�eΪS'�ĶԽ���һ��Ԫ��
	*/
	size_t m = A.msize();
	size_t n = A.nsize();
	int householder_steps = (m == n ? n - 2 : n - 1);
	for (int k = 0; k < householder_steps + 1; ++k) {
		//A��k�д洢uk
		S[k] = 0;
		for (int i = k; i < m; ++i) S[k] = hypot(S[k], A[i][k]);
		if (fabs(S[k]) >0) {
			if (A[k][k] < 0) S[k] = -S[k];
			for (int i = k; i < m; ++i) A[i][k] /= S[k];
			A[k][k] += 1;
		}
		S[k] = -S[k];
		//Hk��A��˴���A��k��֮��
		for (int j = k + 1; j < n; ++j) {
			if (fabs(S[k]) >0) {
				double factor = 0.0;
				for (int i = k; i < m; ++i)factor += A[i][k] * A[i][j];
				factor = -factor / A[k][k];
				for (int i = k; i < m; ++i) { 
					A[i][j] += factor*A[i][k]; 
				}
			}
			e[j] = A[k][j];//e�洢A�Խ�����һ��Ԫ��
		}
		for (int i = k; i < m; ++i)U[i][k] = A[i][k];
		vector<double> work(m);
		if (k < n - 2) {
			e[k] = 0;
			for (int i = k + 1; i < n; ++i)e[k] = hypot(e[k], e[i]);
			if (e[k] != 0)
			{
				if (e[k + 1] < 0)
					e[k] = -e[k];
				for (int i = k + 1; i < n; ++i)
					e[i] /= e[k];
				e[k + 1] += 1;
			}
			e[k] = -e[k];

			if (e[k] != 0) {
				for (int i = k + 1; i < m; ++i)
					work[i] = 0;

				for (int j = k + 1; j < n; ++j)
					for (int i = k + 1; i < m; ++i)
						work[i] += e[j] * A[i][j];

				for (int j = k + 1; j < n; ++j) {
					double t = -e[j] / e[k + 1];
					for (int i = k + 1; i < m; ++i)
						A[i][j] += t * work[i];
				}
			}
			for (int i = k + 1; i < n; ++i)V[i][k] = e[i];
		}
	}

	if (m == n) S[n - 1] = A[n - 1][n - 1];
	if (n > 1)e[n - 2] = A[n - 2][n - 1];
	e[n - 1] = 0;

	//����U
	if (m == n) {
		for (int i = 0; i < m; ++i)
			U[i][n - 1] = 0;
		U[n - 1][n - 1] = 1;
	}
	for (int k = householder_steps; k >= 0; --k) {
		if (fabs(S[k]) > 0) {
			for (int j = k + 1; j < n; ++j) {
				double t = 0;
				for (int i = k; i < m; ++i)t += U[i][k] * U[i][j];
				t = -t / U[k][k];
				for (int i = k; i < m; ++i) {
					U[i][j] += t * U[i][k]; 
				}
			}
			for (int i = k; i < m; ++i)U[i][k] = -U[i][k];
			U[k][k] = 1 + U[k][k];
			for (int i = 0; i < k - 1; ++i)
				U[i][k] = 0;
		}
		else {
			for (int i = 0; i < m; ++i)
				U[i][k] = 0;
			U[k][k] = 1;
		}
	}
	//����V
	for (int k = n - 1; k >= 0; --k) {
		if ((k < n-2) && (e[k] != 0)) {
			for (int j = k + 1; j < n; ++j) {
				double t = 0;
				for (int i = k + 1; i < n; ++i)t += V[i][k] * V[i][j];
				t = -t / V[k + 1][k];
				for (int i = k + 1; i < n; ++i)V[i][j] += t * V[i][k];
			}
		}
		for (int i = 0; i < n; ++i)V[i][k] = 0;
		V[k][k] = 1;
	}
}

int SVD_decomposition(Matrix<double> &A, Matrix<double> &U, vector<double> &S, Matrix<double> &V, double eps){
	/*
		����ֵ�ֽ�
		A=USV��SΪ�ԽǾ���U,VΪ��������
		���أ�A����
	*/
	size_t m = A.msize();
	size_t n = A.nsize();
	vector<double> e(n);
	//householder�ֽ�
	householder_transform(A, U, S, V, e);

	//QR�ֽ�
	int order = n;
	int iter = 0;
	while (order > 0){
		int k;
		int situation;
		
		// situation = 1    e[k+1]-e[order-2]���ܺ��ԣ�e[k]���Ժ��ԣ�s[order - 1]���Ժ���
		// situation = 2     e[k+1]-e[order-2]���ܺ��ԣ�e[k]���Ժ��ԣ�s[k+1]-s[order - 1]�����Ժ��ԣ�s[k]���Ժ���
		// situation = 3     e[k+1]-e[order-2]���ܺ��ԣ�e[k]���Ժ��ԣ�s[k]-s[order - 1]���ܺ���
		// situation = 4     e[order-2]�ɺ���.
		for (k = order - 2; k >= 0; --k){
			if (fabs(e[k]) <= eps*(fabs(S[k]) + fabs(S[k + 1]))){
				//e[k] = 0;
				break;
			}
		}
		if (k == order - 2) situation = 4;
		else{
			int ks;
			for (ks = order - 1; ks > k; --ks){
				double factor = ((ks != order) ? abs(e[ks]) : 0) + ((ks != k + 1) ? abs(e[ks - 1]) : 0);
				if (abs(S[ks]) <= eps*factor){
					//S[ks] = 0;
					break;
				}
			}

			if (ks == k) situation = 3;
			else if (ks == order - 1)situation = 1;
			else{
				situation = 2;
				k = ks;
			}
		}
		k++;
		switch (situation) {
			case 1:
			{
				double f = e[order - 2];
				e[order - 2] = 0;
				for (int j = order - 2; j >= k; --j) {
					double t = hypot(S[j], f);
					double cs = S[j] / t;
					double sn = f / t;
					S[j] = t;
					if (j != k) {
						f = -sn * e[j - 1];
						e[j - 1] = cs * e[j - 1];
					}

					for (int i = 0; i < n; ++i) {
						t = cs*V[i][j] + sn*V[i][order - 1];
						V[i][order - 1] = -sn*V[i][j] + cs*V[i][order - 1];
						V[i][j] = t;
					}
				}
				break;
			}

			case 2:
			{
				double f = e[k - 1];
				e[k - 1] = 0;
				for (int j = k; j < order; ++j){
					double t = hypot(S[j], f);
					double cs = S[j] / t;
					double sn = f / t;
					S[j] = t;
					f = -sn * e[j];
					e[j] = cs * e[j];


					for (int i = 0; i < m; ++i){
						t = cs*U[i][j] + sn*U[i][k - 1];
						U[i][k - 1] = -sn*U[i][j] + cs*U[i][k - 1];
						U[i][j] = t;
					}
				}
				break;
			}

			case 3:
			{
				double scale = max(max(max(max(
					fabs(S[order - 1]), fabs(S[order - 2])), fabs(e[order - 2])),
					fabs(S[k])), fabs(e[k]));
				double sp = S[order - 1] / scale;
				double spm1 = S[order - 2] / scale;
				double epm1 = e[order - 2] / scale;
				double sk = S[k] / scale;
				double ek = e[k] / scale;
				double b = ((spm1 + sp)*(spm1 - sp) + epm1*epm1) / 2.0;
				double c = (sp*epm1) * (sp*epm1);
				double shift = 0;
				if ((b != 0) || (c != 0)){
					shift = sqrt(b*b + c);
					if (b < 0)
						shift = -shift;
					shift = c / (b + shift);
				}
				double f = (sk + sp)*(sk - sp) + shift;
				double g = sk * ek;

				for (int j = k; j < order - 1; ++j){
					double t = hypot(f, g);
					double cs = f / t;
					double sn = g / t;
					if (j != k)e[j - 1] = t;

					f = cs*S[j] + sn*e[j];
					e[j] = cs*e[j] - sn*S[j];
					g = sn * S[j + 1];
					S[j + 1] = cs * S[j + 1];

					for (int i = 0; i < n; ++i)
					{
						t = cs*V[i][j] + sn*V[i][j + 1];
						V[i][j + 1] = -sn*V[i][j] + cs*V[i][j + 1];
						V[i][j] = t;
					}

					t = hypot(f, g);
					cs = f / t;
					sn = g / t;
					S[j] = t;
					f = cs*e[j] + sn*S[j + 1];
					S[j + 1] = -sn*e[j] + cs*S[j + 1];
					g = sn * e[j + 1];
					e[j + 1] = cs * e[j + 1];

					for (int i = 0; i < m; ++i)
					{
						t = cs*U[i][j] + sn*U[i][j + 1];
						U[i][j + 1] = -sn*U[i][j] + cs*U[i][j + 1];
						U[i][j] = t;
					}
				}
				e[order - 2] = f;
				iter = iter + 1;
				break;
			}

			case 4:
			{
				if (S[k] <= 0){
					S[k] = (S[k] < 0) ? -S[k] : 0;
					for (int i = 0; i <= n-1; ++i)
						V[i][k] = -V[i][k];
				}

				while (k < n-1){
					if (S[k] >= S[k + 1])break;

					double t = S[k];
					S[k] = S[k + 1];
					S[k + 1] = t;

					for (int i = 0; i < n; ++i)
						swap(V[i][k], V[i][k + 1]);

					for (int i = 0; i < m; ++i)
						swap(U[i][k], U[i][k + 1]);
					k++;
				}
				iter = 0;
				order--;
				break;
			}
		}
	}
	//������
	size_t N = S.size();
	double factor = N*S[0] * eps;
	int r = 0;
	for (int i = 0; i < N; ++i) {
		if (S[i] > factor) ++r;
		//else S[i] = 0;
	}
	return r;
}

vector<shared_ptr<Linear_solution<double>>> solve_linear_SVD_method(const Matrix<double> &A,const Matrix<double> &b, double eps) {
	/*
		���Է�����ͨ�ýⷨ
		������
		A��m*mϵ������
		b��m*n�Ҷ���
		���أ�
		��
	*/
	size_t m = A.msize();
	size_t n = A.nsize();
	size_t p = max(m, n);
	Matrix<double> U(p, p);
	Matrix<double> V(p, p);
	Matrix<double> B(p,p);
	vector<double> S(p);
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			B[i][j] = A[i][j];
		}
	}
	//��A����ֵ�ֽ�ΪUSV
	int r = SVD_decomposition(B, U, S, V, eps);
	double condition_number = 0;
	if (r == n) condition_number = S[0] / S[n - 1];

	//�������
	vector<shared_ptr<Linear_solution<double>>>solutions;
	size_t question_num = b.nsize();
	size_t general_solution_num = n - r;
	for (int question = 0; question < question_num; ++question) {
		shared_ptr<Linear_solution<double>> solution(new Linear_solution<double>(n, general_solution_num));
		vector<double> temp(n, 0);
		for (int j = 0; j < r; ++j) {
			double factor = 0.0;
			if (fabs(S[j]) > 0) {
				for (int i = 0; i < m; ++i) factor += U[i][j] * b[i][question];
				factor /= S[j];
			}
			temp[j] = factor;
		}
		vector<double> particular_solution(n,0);
		for (int i = 0; i < n; ++i) {
			double factor = 0.0;
			for (int j = 0; j < r; ++j)factor += V[i][j] * temp[j];
			particular_solution[i] = factor;
		}
		Matrix<double> general_solution(general_solution_num, n);
		if (general_solution_num > 0) {
			for (int j = r; j < n; ++j) {
				for (int i = 0; i < n; ++i) {
					general_solution[j - r][i] = V[i][j];
				}
			}
		}
		vector<double> b_solved(m,0);
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n;++j ) {
				b_solved[i] += A[i][j]*particular_solution[j];
			}
		}
		bool solution_exist = true;
		for (int i = 0; i < m; ++i) {
			if (fabs(b[i][question] - b_solved[i]) > EPS) {
				solution_exist = false;
				break;
			}
		}
		if (solution_exist) {
			solution->set_solution_value(particular_solution, general_solution);
		}
		else {
			solution->set_solution_value(particular_solution, Matrix<double>());
			solution->set_solution_type(!solution_exist);
		}
		solutions.push_back(solution);
	}

	return solutions;
}

shared_ptr<vector<double>> solve_tridiagonal(const vector<double> &a,const vector<double> &b,const vector<double> &c,const vector<double> &r) {
	/*
		�����Խ�ϵ�����󷽳�
		������
		a���Խ�����һ��Ԫ�أ�n-1ά����
		b���Խ���Ԫ�أ�nά����
		c���Խ�����һ��Ԫ�أ�n-1ά����
		���أ�
		��
	*/
	size_t n = b.size();
	shared_ptr<vector<double>> solution(new vector<double>(n));
	shared_ptr<vector<double>> solution_null;
	if (b.size() != a.size() + 1 || a.size() != c.size()) {
		cout << "���Խ������ʽ����" << endl;
		return solution_null;
	}
	else if (fabs(b[0]) < EPS) {
		cout << "b[0]����Ϊ0" << endl;
		return solution_null;
	}
	vector<double> temp(n,0);
	double factor = b[0];
	(*solution)[0] = r[0] / factor;
	for (int i = 1; i < n; ++i) {
		temp[i] = c[i-1]/factor;
		factor = b[i] - a[i-1] * temp[i];
		if (fabs(factor) < EPS) {
			cout << "��ֵ���ȶ����޷�����" << endl;
			return solution_null;
		}
		(*solution)[i] = (r[i] - a[i-1] * (*solution)[i - 1]) / factor;
	}
	for (int i = n - 2; i >= 0; --i) {
		(*solution)[i] -= temp[i + 1] * (*solution)[i + 1];
	}
	return solution;
}

shared_ptr<vector<double>> solve_cyclic_tridiagonal(const vector<double> &a, const vector<double> &b, const vector<double> &c,const double alpha,const double beta, const vector<double> &r) {
	/*
		���������Խ�ϵ�����󷽳�Ax=r
		������
		a���Խ�����һ��Ԫ�أ�n-1ά����
		b���Խ���Ԫ�أ�nά����
		c���Խ�����һ��Ԫ�أ�n-1ά����
		alpha�����½�Ԫ��
		beta�����Ͻ�Ԫ��
		���أ�
		��
	*/
	size_t n = b.size();
	shared_ptr<vector<double>> solution_null;
	if (b.size() != a.size() + 1 || a.size() != c.size()) {
		cout << "���Խ������ʽ����" << endl;
		return solution_null;
	}
	else if (n<=2) {
		cout << "�����С���޷�����" << endl;
		return solution_null;
	}
	double factor, gamma;
	vector<double> new_b(n), u(n,0), z(n);
	gamma = -b[0];
	new_b[0] = b[0] - gamma;
	new_b[n - 1] = b[n - 1] - alpha*beta / gamma;
	for (int i = 1; i < n - 1; ++i) new_b[i] = b[i];
	auto solution_x = solve_tridiagonal(a, new_b, c, r);//���Ay=r
	if (solution_x.get() == nullptr) return solution_x;
	u[0] = gamma;
	u[n - 1] = alpha;
	auto solution_z = solve_tridiagonal(a, new_b, c, u);//���Az=u
	if (solution_z.get() == nullptr) return solution_z;
	factor = ((*solution_x)[0] + beta*(*solution_x)[n - 1] / gamma) / (1.0 + (*solution_z)[0] + beta*(*solution_z)[n - 1] / gamma);
	for (int i = 0; i < n; ++i) (*solution_x)[i] -= factor*(*solution_z)[i];//���x
	return solution_x;
}

shared_ptr<vector<double>> solve_Vandermonde(const vector<double> &A, const vector<double> &b) {
	/*
		��Vandermondeϵ�����󷽳���Ax=b
		ϵ���������磺
		1	x0			x0^2		...	x0^(n-1)
		1	x1			x1^2		...	x1^(n-1)
		...
		1	x(n-1)	x(n-1)^2	...	x(n-1)^(n-1)

		������
		A��[x0	x1	...	x(n-1)]����
		b���Ҷ�����
		���أ�
		��
	*/
	size_t n = A.size();
	shared_ptr<vector<double>> solution;
	if (n != b.size()) {
		cout << "��ʽ����" << endl;
		return solution;
	}
	solution.reset(new vector<double>(n,0));
	vector<double> temp(n, 0);
	temp[n - 1] = -A[0];
	for (int i = 1; i < n; ++i) {
		for (int j = (n - 1 - i); j < n - 1; ++j)temp[j] -= A[i]*temp[j + 1];
		temp[n - 1] -= A[i];
	}
	for (int i = 0; i < n; ++i) {
		double factor = (double)n;
		for (int j = n - 1; j > 0; --j) factor = j*temp[j] + A[i] * factor;
		factor = b[i] / factor;
		double r = 1.0;
		for (int j = n - 1; j >= 0; --j) {
			(*solution)[j] += r*factor;
			r = temp[j] + A[i] * r;
		}
	}
	return solution;
}

shared_ptr<vector<double>> solve_Vandermonde_transposition(const vector<double> &A, const vector<double> &b) {
	/*
		��Vandermondeϵ�����󷽳���Ax=b
		ϵ���������磺
		1				1				...		1
		x0				x1				...		x(n-1)
		x0^2		x1^2		...		x(n-1)^2
		...
		x0^(n-1)	x1^(n-1)	...		x(n-1)^(n-1)

		������
		A��[x0	x1	...	x(n-1)]����
		b���Ҷ�����
		���أ�
		��
	*/
	size_t n = A.size();
	shared_ptr<vector<double>> solution;
	if (n != b.size()) {
		cout << "��ʽ����" << endl;
		return solution;
	}
	solution.reset(new vector<double>(n));
	vector<double> temp(n,0);
	double factor;
	if (n == 1) (*solution)[0] = b[0];
	else {
		temp[n - 1] = -A[0];
		for (int i = 1; i < n; ++i) {
			factor = -A[i];
			for (int j = (n - 1 - i); j < n - 1; ++j)temp[j] += factor*temp[j + 1];
			temp[n - 1] += factor;
		}

		for (int i = 0; i < n; ++i) {
			factor = A[i];
			double denominator = 1.0, r = 1.0;
			double numerator = b[n - 1];
			for (int j = n - 1; j > 0; --j) {
				r = temp[j] + factor*r;
				numerator += b[j - 1] * r;
				denominator = factor*denominator + r;
			}
			(*solution)[i] = numerator / denominator;
		}
	}
	return solution;
}