#pragma once
#include <vector>
#include<iostream>
#include<memory>
#include<algorithm>
#include"Matrix.hpp"

template <typename value_type>
class Linear_solution
{
private:
	std::vector<value_type> *particular_solution=nullptr;
	std::vector<std::vector<value_type>> *general_solution = nullptr;//每个解按行存储
	size_t solution_dimension = 0;
	size_t general_solution_num = 0;
	bool least_square_solution = false;
public:
	//构造函数
	Linear_solution();
	Linear_solution(size_t solution_dimension, size_t general_solution_num);
	Linear_solution(const Linear_solution<value_type> &input);//Linear_solution复制构造
	Linear_solution(const std::vector<value_type> &particular_solution, const Matrix<value_type> &general_solution);//特解和通解构造

	//赋值
	Linear_solution<value_type> & operator=(const Linear_solution<value_type> &input);//=Linear_solution

	//读取成员
	inline int get_solution_dimension() const;
	inline int get_general_solution_num() const;
	inline bool is_least_square_solution() const;
	inline const std::vector<value_type> * get_particular_solution() const;
	inline const std::vector<std::vector<value_type>> * get_general_solution() const;

	//设置成员
	inline void set_solution_value(const std::vector<value_type> &particular_solution, const Matrix<value_type> &general_solution);
	inline void set_solution_type(bool type);
	//打印
	inline void print() const;

	//析构
	~Linear_solution();
};

template <typename value_type>
Linear_solution<value_type>::Linear_solution() {
	particular_solution = nullptr;
	general_solution = nullptr;
}

template <typename value_type>
Linear_solution<value_type>::Linear_solution(size_t solution_dimension, size_t general_solution_num) {
	this->solution_dimension = solution_dimension;
	this->general_solution_num = general_solution_num;
	particular_solution = new std::vector<value_type>(solution_dimension,0);
	general_solution = new std::vector<std::vector<value_type>>();
	for (int i = 0; i < general_solution_num; ++i) {
		std::vector<value_type> line(solution_dimension, 0);
		(*general_solution).push_back(line);
	}
}

template <typename value_type>
Linear_solution<value_type>::Linear_solution(const Linear_solution<value_type> &input) {
	if (this->solution_dimension != input.get_solution_dimension() || this->general_solution_num != input.get_general_solution_num()) {
		if (particular_solution != nullptr) delete particular_solution;
		if (general_solution != nullptr)delete general_solution;
		this->solution_dimension = input.get_solution_dimension();
		this->general_solution_num = input.get_general_solution_num();
		particular_solution = new std::vector<value_type>(solution_dimension, 0);
		general_solution = new std::vector<std::vector<value_type>>();
		for (int i = 0; i < general_solution_num; ++i) {
			std::vector<value_type> line(solution_dimension, 0);
			(*general_solution).push_back(line);
		}
	}
	const std::vector<value_type> *input_particular_solution = input.get_particular_solution();
	const std::vector<std::vector<value_type>> *input_general_solution = input.get_general_solution();
	for (int i = 0; i < solution_dimension; ++i) {
		(*particular_solution)[i] = (*input_particular_solution)[i];
	}
	for (int i = 0; i < general_solution_num; ++i) {
		for (int j = 0; j < solution_dimension; ++j) {
			(*general_solution)[i][j] = (*input_general_solution)[i][j];
		}
	}
}

template <typename value_type>
Linear_solution<value_type>::Linear_solution(const std::vector<value_type> &particular_solution, const Matrix<value_type> &general_solution) {
	int solution_dimension = particular_solution.size();
	int general_solution_num = general_solution.msize();
	if (general_solution_num > 0 && solution_dimension != general_solution.nsize()) {
		cout << "输入格式错误" << endl;
		return;
	}
	if (this->solution_dimension != solution_dimension || this->general_solution_num != general_solution_num) {
		if (this->particular_solution != nullptr) delete this->particular_solution;
		if (this->general_solution != nullptr)delete this->general_solution;
		this->solution_dimension = solution_dimension;
		this->general_solution_num = general_solution_num;
		this->particular_solution = new std::vector<value_type>(solution_dimension, 0);
		this->general_solution = new std::vector<std::vector<value_type>>();
		for (int i = 0; i < general_solution_num; ++i) {
			std::vector<value_type> line(solution_dimension, 0);
			(*this->general_solution).push_back(line);
		}
	}
	for (int i = 0; i < solution_dimension; ++i) {
		(*this->particular_solution)[i] = particular_solution[i];
	}
	for (int i = 0; i < general_solution_num; ++i) {
		for (int j = 0; j < solution_dimension; ++j) {
			(*this->general_solution)[i][j] = general_solution[i][j];
		}
	}
}

template <typename value_type>
Linear_solution<value_type> & Linear_solution<value_type>::operator=(const Linear_solution<value_type> &input) {
	if (this->solution_dimension != input.get_solution_dimension() || this->general_solution_num != input.get_general_solution_num()) {
		if (particular_solution != nullptr) delete particular_solution;
		if (general_solution != nullptr)delete general_solution;
		this->solution_dimension = input.get_solution_dimension();
		this->general_solution_num = input.get_general_solution_num();
		particular_solution = new std::vector<value_type>(solution_dimension, 0);
		general_solution = new std::vector<std::vector<value_type>>();
		for (int i = 0; i < general_solution_num; ++i) {
			std::vector<value_type> line(solution_dimension, 0);
			(*general_solution).push_back(line);
		}
	}
	const std::vector<value_type> *input_particular_solution = input.get_particular_solution();
	const std::vector<std::vector<value_type>> *input_general_solution = input.get_general_solution();
	for (int i = 0; i < solution_dimension; ++i) {
		(*particular_solution)[i] = (*input_particular_solution)[i];
	}
	for (int i = 0; i < general_solution_num; ++i) {
		for (int j = 0; j < solution_dimension; ++j) {
			(*general_solution)[i][j] = (*input_general_solution)[i][j];
		}
	}
	return *this;
}

template <typename value_type>
inline int Linear_solution<value_type>::get_solution_dimension() const {
	return solution_dimension;
}

template <typename value_type>
inline int Linear_solution<value_type>::get_general_solution_num() const {
	return general_solution_num;
}

template <typename value_type>
inline bool Linear_solution<value_type>::is_least_square_solution() const {
	return least_square_solution;
}

template <typename value_type>
inline const std::vector<value_type> * Linear_solution<value_type>::get_particular_solution() const {
	return particular_solution;
}

template <typename value_type>
inline const std::vector<std::vector<value_type>> * Linear_solution<value_type>::get_general_solution() const {
	return general_solution;
}

template <typename value_type>
inline void Linear_solution<value_type>::set_solution_value(const std::vector<value_type> &particular_solution,const Matrix<value_type> &general_solution) {
	size_t solution_dimension = particular_solution.size();
	size_t general_solution_num = general_solution.msize();
	if (general_solution_num > 0 && solution_dimension != general_solution.nsize()) {
		cout << "输入格式错误" << endl;
		return;
	}
	if (this->solution_dimension != solution_dimension || this->general_solution_num != general_solution_num) {
		if (this->particular_solution != nullptr) delete this->particular_solution;
		if (this->general_solution != nullptr)delete this->general_solution;
		this->solution_dimension = solution_dimension;
		this->general_solution_num = general_solution_num;
		this->particular_solution = new std::vector<value_type>(solution_dimension, 0);
		this->general_solution = new std::vector<std::vector<value_type>>();
		for (int i = 0; i < general_solution_num; ++i) {
			std::vector<value_type> line(solution_dimension, 0);
			(*this->general_solution).push_back(line);
		}
	}
	for (int i = 0; i < solution_dimension; ++i) {
		(*this->particular_solution)[i] = particular_solution[i];
	}
	for (int i = 0; i < general_solution_num; ++i) {
		for (int j = 0; j < solution_dimension; ++j) {
			(*this->general_solution)[i][j] = general_solution[i][j];
		}
	}
}

template <typename value_type>
inline void Linear_solution<value_type>::set_solution_type(bool type) {
	least_square_solution = type;
}

template <typename value_type>
inline void Linear_solution<value_type>::print() const {
	if (least_square_solution) {
		std::cout << "无解，给出最小二乘解：" << std::endl;
		std::cout << "[ ";
		if (particular_solution != nullptr) {
			for (auto &elem : *particular_solution) {
				std::cout << elem << " ";
			}
		}
		std::cout << "]" << std::endl;
		return;
	}
	if (general_solution == nullptr || general_solution->size() == 0) {
		std::cout << "解为：" << std::endl << "[ ";
	}
	else {
		std::cout << "特解为：" << std::endl << "[ ";
	}
	if (particular_solution != nullptr) {
		for (auto &elem : *particular_solution) {
			std::cout << elem << " ";
		}
	}
	std::cout << "]" << std::endl;
	if (general_solution != nullptr&&general_solution->size()!=0) {
		std::cout << "通解为：" << std::endl;
		for (auto &solution : *general_solution) {
			std::cout << "[ ";
			for (auto &elem : solution) {
				std::cout << elem << " ";
			}
			std::cout << "]" << std::endl;
		}
	}
}

template <typename value_type>
Linear_solution<value_type>::~Linear_solution() {
	if (particular_solution != nullptr) delete particular_solution;
	if (general_solution != nullptr)delete general_solution;
}
