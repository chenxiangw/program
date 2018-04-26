#include<iostream>
#include<vector>
#include<string>
#include<map>
#include <memory>  
#include <functional> 
#include"interaction.h"
#include"Matrix.h"
#include"numerical_recipes.h"
using namespace std;

map < string, shared_ptr<Matrix<double>>> matrixes;
map < string, unique_ptr<Function> > commands;
map < string, string > function_instruction;

template<typename ...Arguments>
void call(const string &input, Arguments&& ... args){
	const Function &base = *commands[input];
	const function<void(Arguments...)> & fun = static_cast<const Function_type<Arguments...>&>(base).callback;
	fun(forward<Arguments>(args)...);
}

void interaction() {

	unique_ptr<Function_type<const string>> func_dict(new Function_type<const string>(&dict));
	unique_ptr<Function_type<const string>> func_matrix(new Function_type<const string>(&matrix));
	unique_ptr<Function_type<const string>> func_inverse(new Function_type<const string>(&inverse));
	unique_ptr<Function_type<const string>> func_determinant(new Function_type<const string>(&determinant));
	unique_ptr<Function_type<const string>> func_solve(new Function_type<const string>(&solve));
	commands.insert(map < string, unique_ptr<Function> >::value_type("dict", move(func_dict)));
	commands.insert(map < string, unique_ptr<Function> >::value_type("matrix", move(func_matrix)));
	commands.insert(map < string, unique_ptr<Function> >::value_type("inverse", move(func_inverse)));
	commands.insert(map < string, unique_ptr<Function> >::value_type("determinant", move(func_determinant)));
	commands.insert(map < string, unique_ptr<Function> >::value_type("solve", move(func_solve)));
	function_instruction.insert(map < string, string >::value_type("dict", "dict(查询所有指令)\ndict 指令名(查询该指令)\n"));
	function_instruction.insert(map < string, string >::value_type("matrix", "matrix(查询已存在矩阵)\nmatrix 矩阵名(查询或创建该矩阵)\n"));
	function_instruction.insert(map < string, string >::value_type("inverse", "inverse 矩阵名(求矩阵的逆)\n"));
	function_instruction.insert(map < string, string >::value_type("determinant", "determinant 矩阵名(求矩阵的行列式)\n"));
	function_instruction.insert(map < string, string >::value_type("solve", "solve 系数矩阵A 列向量矩阵b(求方程的解)\n"));

	cout << "请输入指令（dict查询指令；quit退出）：" << endl;
	string instruction="";
	while (true) {
		getline(cin, instruction);
		if (instruction == "quit") break;
		size_t sep = instruction.find_first_of(' ');
		string command = instruction.substr(0, sep);
		string parameter = "";
		if(sep!=string::npos) parameter=instruction.substr(sep + 1, instruction.size());
		if (commands.find(command) != commands.end()) {
			call(command, parameter);
			cout << "请输入指令（dict查询指令；quit退出）：" << endl;
		}
		else {
			cout << "指令无法识别，请重新输入：" << endl;
		}

	}
}