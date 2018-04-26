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
	function_instruction.insert(map < string, string >::value_type("dict", "dict(��ѯ����ָ��)\ndict ָ����(��ѯ��ָ��)\n"));
	function_instruction.insert(map < string, string >::value_type("matrix", "matrix(��ѯ�Ѵ��ھ���)\nmatrix ������(��ѯ�򴴽��þ���)\n"));
	function_instruction.insert(map < string, string >::value_type("inverse", "inverse ������(��������)\n"));
	function_instruction.insert(map < string, string >::value_type("determinant", "determinant ������(����������ʽ)\n"));
	function_instruction.insert(map < string, string >::value_type("solve", "solve ϵ������A ����������b(�󷽳̵Ľ�)\n"));

	cout << "������ָ�dict��ѯָ�quit�˳�����" << endl;
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
			cout << "������ָ�dict��ѯָ�quit�˳�����" << endl;
		}
		else {
			cout << "ָ���޷�ʶ�����������룺" << endl;
		}

	}
}