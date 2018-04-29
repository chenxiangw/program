#include<iostream>
#include<string>
#include<map>
#include<vector>
#include <functional> 
#include"Matrix.h"
#include"interaction.h"
#include"numerical_recipes.h"
using namespace std;
extern map < string, shared_ptr<Matrix<double>> > matrixes;
extern map < string, unique_ptr<Function> > commands;
extern map < string, string > function_instruction;

bool split_to_double(const string &line, vector<double> &nums, char sep = ' ') {
	string temp = "";
	int size = line.size();
	for (int pos = 0; pos < size; ++pos) {
		if (line[pos] == sep&&temp.size() != 0) {
			double num = stod(temp);
			nums.push_back(num);
			temp = "";
		}
		else {
			temp += line[pos];
		}
	}
	if (temp.size() != 0) nums.push_back(stod(temp));
	return true;
}

void split_to_string(const string &line, vector<string> &nums, char sep = ' ') {
	string temp = "";
	int size = line.size();
	for (int pos = 0; pos < size; ++pos) {
		if (line[pos] == sep&&temp.size() != 0) {
			nums.push_back(temp);
			temp = "";
		}
		else {
			temp += line[pos];
		}
	}
	if (temp.size() != 0) nums.push_back(temp);
}

void dict(const string &input) {
	if (input.size() == 0) {
		cout << "指令操作有：" << endl;
		for (auto &elem : commands) {
			cout << elem.first << " ";
		}
		cout << endl;
	}
	else {
		if (function_instruction.find(input) != function_instruction.end()) {
			cout << input << "指令说明：" << endl;
			cout << function_instruction[input];
		}
		else {
			cout << "没有找到该指令" << endl;
		}
	}
}

void matrix(const string &input) {
	if (input.size() == 0) {
		cout << "当前矩阵有：" << endl;
		for (auto &elem : matrixes) {
			cout << elem.first << " ";
		}
		cout << endl;
	}
	else if(input.find_first_of(' ')!=string::npos){
		cout << "矩阵名称错误：不能包含空格" << endl;
	}
	else if (matrixes.find(input) != matrixes.end()) {
		cout << "矩阵" << input << "为：" << endl;
		matrixes[input]->print();
	}
	else {
		cout << "创建矩阵" << input << "（空格分隔，输入quit放弃，ok完成）：" << endl;
		vector<vector<double>> temp;
		shared_ptr<Matrix<double>> new_matrix(new Matrix<double>());
		int n = -1;
		string line = "";
		while (true) {
			line = "";
			getline(cin, line);
			if (line == "quit") break;
			if (line == "ok") {
				*new_matrix = temp;
				matrixes.insert(map < string, shared_ptr<Matrix<double>> >::value_type(input, new_matrix));
				break;
			}
			else {
				vector<double> nums;
				bool succeed = split_to_double(line, nums);
				if (!succeed) {
					cout << "输入错误，请重新输入：" << endl;
					continue;
				}
				if (n == -1) n = nums.size();
				else {
					if (n != nums.size()) {
						cout << "输入错误，请重新输入：" << endl;
						continue;
					}
				}
				temp.push_back(nums);
			}
		}
	}
}

void clear(const string &input) {
	if (input.size() == 0) {
		cout << "是否清除所有矩阵？(Y/N)" << endl;
		string confirm="";
		while (true) {
			confirm = "";
			getline(cin, confirm);
			if (confirm.size() > 0 && toupper(confirm[0]) == 'Y') {
				matrixes.clear();
				break;
			}
			else if (confirm.size() > 0 && toupper(confirm[0]) == 'N') {
				break;
			}
			else {
				cout << "输入错误，请重新输入：" << endl;
				continue;
			}
		}
	}
	else if (matrixes.find(input) != matrixes.end()) {
		matrixes.erase(input);
	}
	else {
		cout << "删除错误：未找到该矩阵" << endl;
	}
}

void inverse(const std::string &input){
	if (matrixes.find(input) != matrixes.end()) {
		auto result = matrixes[input]->inverse();
		if (result.get() != nullptr){
			cout << "逆矩阵为：" << endl;
			result->print();
		}
	}
	else {
		cout << "参数错误" << endl;
	}
}

void determinant(const std::string &input){
	if (matrixes.find(input) != matrixes.end()) {
		cout << "行列式为：";
		cout<<matrixes[input]->determinant()<<endl;
	}
	else {
		cout << "参数错误" << endl;
	}
}

void solve(const std::string &input){
	vector<string> parameters;
	split_to_string(input, parameters);
	if (parameters.size() != 2) {
		cout << "参数错误" << endl;
		return;
	}
	for (auto &elem : parameters) {
		if (matrixes.find(elem) == matrixes.end()) {
			cout << "矩阵不存在" << endl;
			return;
		}
	}
	if (matrixes[parameters[0]]->msize() != matrixes[parameters[1]]->msize()) {
		cout << "矩阵格式错误" << endl;
		return;
	}
	auto result = solve_linear_SVD_method(*matrixes[parameters[0]], *matrixes[parameters[1]]);
	int question_num = result.size();
	for (int i = 0; i < question_num; ++i) {
		cout << "第" << i+1 << "个方程：" << endl;
		result[i]->print();
	}
}

