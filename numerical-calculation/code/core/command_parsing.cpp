#include<iostream>
#include<string>
#include<map>
#include<utility>
#include<vector>
#include<limits>
#include <stdexcept>  
#include <boost\spirit\include\qi.hpp>
#include<boost\spirit\include\phoenix.hpp>
#include<boost\bind\bind.hpp>
#include <boost\lambda\lambda.hpp>
#include"Matrix.h"
#include"Evaluation_function.h"
#include"numerical_recipes.h"
#include"builtin_functions.h"


std::map < std::string, std::string > command_instruction;
std::map < std::string, std::shared_ptr<double>> symbols;
std::map < std::string, std::string> functions;
std::map < std::string, std::shared_ptr<Matrix<double>>> matrixes;
typedef  std::pair<void(*)(const std::vector<std::string> *), int> function_type;
std::map < std::string, function_type > builtin_functions;
std::map < std::string, void(*)(const std::vector<std::string> *input)> commands;

std::vector<std::string> command_params;
std::string command_name;
std::vector<std::string> function_params;
std::vector<std::string> function_names;
std::string var_name;
std::vector<std::shared_ptr<Matrix<double>>> result_matrix;

bool split_to_double(const std::string &line, std::vector<double> &nums, char sep = ' ') {
	std::string temp = "";
	int size = line.size();
	try {
		for (int pos = 0; pos < size; ++pos) {
			if (line[pos] == sep && temp.size() != 0) {
				double num;
				num = stod(temp);
				nums.push_back(num);
				temp = "";
			}
			else {
				temp += line[pos];
			}
		}
		if (temp.size() != 0) nums.push_back(stod(temp));
	}
	catch (const std::invalid_argument& ia) {
		return false;
	}
	return true;
}

void dict(const std::vector<std::string> *input) {
	int n = input->size();
	if (n == 0) {
		std::cout << "指令操作有：" << std::endl;
		for (auto &elem : command_instruction) {
			std::cout << elem.first << " ";
		}
		std::cout << std::endl;
	}
	else {
		for (int i = 0; i < n; ++i) {
			if (command_instruction.find((*input)[i]) != command_instruction.end()) {
				std::cout << (*input)[i] << "指令说明：" << std::endl;
				std::cout << command_instruction[(*input)[i]];
			}
			else {
				std::cout << "没有找到指令" << (*input)[i] << std::endl;
			}
		}
	}
}

void var(const std::vector<std::string> *input) {
	using namespace std;
	int n = input->size();
	if (n == 0) {
		cout << "当前变量有：" << endl;
		for (auto &elem : symbols) {
			cout << elem.first << " ";
		}
		cout << endl;
	}
	else {
		for (int i = 0; i < n; ++i) {
			if (symbols.find((*input)[i]) == symbols.end()) {
				if (matrixes.find((*input)[i]) != matrixes.end()) {
					cout << "已存在同名矩阵，请重新命名" << endl;
					continue;
				}
				if (builtin_functions.find((*input)[i]) != builtin_functions.end()) {
					cout << "与内建函数名冲突" << endl;
					continue;
				}
				if (functions.find((*input)[i]) != functions.end()) {
					cout << "已存在同名函数，请重新命名" << endl;
					continue;
				}
				symbols.insert(map < string, shared_ptr<double>>::value_type((*input)[i], nullptr));
			}
			else {
				if (symbols[(*input)[i]] != nullptr) {
					cout << (*input)[i] << " = " << *symbols[(*input)[i]] << endl;
				}
				else {
					cout << (*input)[i] << endl;
				}
			}
		}
	}
}

void func(const std::vector<std::string> *input) {
	using namespace std;
	int n = input->size();
	if (n == 0) {
		cout << "当前函数有：" << endl;
		for (auto &elem : functions) {
			cout << elem.first << " ";
		}
		cout << endl;
	}
	else {
		for (int i = 0; i < n; ++i) {
			if (functions.find((*input)[i]) != functions.end()) {
				cout << "函数" << (*input)[i] << "为：" << endl;
				cout << functions[(*input)[i]] << endl;
			}
			else if (matrixes.find((*input)[i]) != matrixes.end()) {
				cout << "已存在同名矩阵，请重新命名"<< endl;
				break;
			}
			else if (builtin_functions.find((*input)[i]) != builtin_functions.end()) {
				cout << "与内建函数名冲突" << endl;
				break;
			}
			else if (symbols.find((*input)[i]) != symbols.end()) {
				cout << "已存在同名变量，请重新命名" << endl;
				break;
			}
			else{
				cout << "创建函数" << (*input)[i] << "（变量名为x，输入quit放弃，ok完成）：" << endl;
				string new_function = "";
				int n = -1;
				string line = "";
				while (true) {
					line = "";
					getline(cin, line);
					if (line == "quit") break;
					if (line == "ok") {
						functions.insert(map < string, string >::value_type((*input)[i], new_function));
						break;
					}
					else {
						new_function += line;
					}
				}
			}
		}
	}
}

void matrix(const std::vector<std::string> *input) {
	using namespace std;
	int num = input->size();
	if (num == 0) {
		cout << "当前矩阵有：" << endl;
		for (auto &elem : matrixes) {
			cout << elem.first << " ";
		}
		cout << endl;
	}
	else {
		for (int i = 0; i < num; ++i) {
			if (matrixes.find((*input)[i]) != matrixes.end()) {
				cout << "矩阵" << (*input)[i] << "为：" << endl;
				matrixes[(*input)[i]]->print();
			}
			else if (builtin_functions.find((*input)[i]) != builtin_functions.end()) {
				cout << "与内建函数名冲突" << endl;
				break;
			}
			else if (functions.find((*input)[i]) != functions.end()) {
				cout << "已存在同名函数，请重新命名" << endl;
				break;
			}
			else if (symbols.find((*input)[i]) != symbols.end()) {
				cout << "已存在同名变量，请重新命名" << endl;
				break;
			}
			else {
				cout << "创建矩阵" << (*input)[i] << "（空格分隔，输入quit放弃，ok完成）：" << endl;
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
						matrixes.insert(map < string, shared_ptr<Matrix<double>> >::value_type((*input)[i], new_matrix));
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
	}
}

void clear(const std::vector<std::string> *input) {
	using namespace std;
	int n = input->size();
	if (n == 0) {
		cout << "是否清空所有变量?(Y/N)" << endl;
		string confirm;
		while (true) {
			confirm = "";
			getline(cin, confirm);
			if (confirm.size() > 0 && toupper(confirm[0]) == 'Y') {
				matrixes.clear();
				functions.clear();
				symbols.clear();
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
	else {
		for (int i = 0; i < n; ++i) {
			if (functions.find((*input)[i]) != functions.end()) {
				functions.erase((*input)[i]);
			}
			else if (matrixes.find((*input)[i]) != matrixes.end()) {
				matrixes.erase((*input)[i]);
			}
			else if (symbols.find((*input)[i]) != symbols.end()) {
				symbols.erase((*input)[i]);
			}
			else {
				cout << "没有找到变量" << (*input)[i] << endl;
			}
		}
	}

}

bool run_command(const std::string &command_name,const std::vector<std::string> *params) {
	using namespace std;
	if (commands.find(command_name) != commands.end()) {
		commands[command_name](params);
	}
	else if(symbols.find(command_name)!= symbols.end()){
		cout << command_name << " = " << *symbols[command_name] << endl;
	}
	else if(functions.find(command_name) != functions.end()){
		cout << "函数" << command_name << "为：" << endl;
		cout << functions[command_name] << endl;
	}
	else if (matrixes.find(command_name) != matrixes.end()) {
		cout << "矩阵" << command_name << "为：" << endl;
		matrixes[command_name]->print();
	}
	else {
		return false;
	}
	return true;
}

void run_function(std::vector<std::string> *function_names,std::vector<std::string> *params) {
	using namespace std;
	if (function_names->size() < 1) {
		throw  string("指令错误") ;
	}
	string function_name = function_names->back();
	function_names->pop_back();
	if (builtin_functions.find(function_name) != builtin_functions.end()) {
		int n = builtin_functions[function_name].second;
		if (n >= params->size() || (*params)[params->size() - n - 1] != "(") {
			throw  string("参数错误");
		}
		vector<string> now_params(n);
		for (int i = n - 1; i >= 0; --i) {
			if (symbols.find(params->back()) != symbols.end()) {
				now_params[i] = to_string(*symbols[params->back()]);
			}
			else {
				now_params[i] = params->back();
			}
			params->pop_back();
		}
		params->pop_back();
		builtin_functions[function_name].first(&now_params);
	}
	else if(functions.find(function_name) != functions.end()){
		Evaluation_function func(function_name);
		if (params->size() <= 1 || (*params)[params->size() - 2] != "(") {
			throw   string("参数错误");
		}
		double x;
		if (symbols.find(params->back()) != symbols.end()) {
			x = *symbols[params->back()];
		}
		else {
			x = stod(params->back());
		}
		params->pop_back();
		params->pop_back();
		function_params.push_back(to_string(func(x)));
	}
	else {
		string error = function_name + "不存在";
		throw error;
	}

}

void init_commands(){
	using namespace std;
	command_instruction.insert(map < string, string >::value_type("dict", "dict：查询所有指令\ndict 指令名：查询该指令\n"));
	commands.insert(map < string, void(*)(const vector<string> *input)>::value_type("dict", dict));

	command_instruction.insert(map < string, string >::value_type("var", "var：查询已存在变量\nvar 变量名：查询或创建变量\n"));
	commands.insert(map < string, void(*)(const vector<string> *input)>::value_type("var", var));

	command_instruction.insert(map < string, string >::value_type("matrix", "matrix：查询已存在矩阵\nmatrix 矩阵名：查询或创建矩阵\n"));
	commands.insert(map < string, void(*)(const vector<string> *input)>::value_type("matrix", matrix));

	command_instruction.insert(map < string, string >::value_type("function", "function：查询已存在函数\nclear 函数名：查询或创建函数\n"));
	commands.insert(map < string, void(*)(const vector<string> *input)>::value_type("function", func));

	command_instruction.insert(map < string, string >::value_type("clear", "clear：删除所有矩阵/函数\nclear 矩阵/函数名：删除该矩阵/函数\n"));
	commands.insert(map < string, void(*)(const vector<string> *input)>::value_type("clear", clear));

	command_instruction.insert(map < string, string >::value_type("inverse", "inverse(矩阵名)：求矩阵的逆\n"));
	builtin_functions.insert(pair<string, function_type>("inverse", function_type(inverse, 1)));

	command_instruction.insert(map < string, string >::value_type("determinant", "determinant(矩阵名)：求矩阵的行列式\n"));
	builtin_functions.insert(pair<string, function_type>("determinant", function_type(determinant, 1)));

	command_instruction.insert(map < string, string >::value_type("solve", "solve(系数矩阵A,列向量矩阵b)：求方程的解\n"));
	builtin_functions.insert(pair<string, function_type>("solve", function_type(solve, 2)));

	command_instruction.insert(map < string, string >::value_type("integrate", "integrate(函数名, 积分下界, 积分上界)：求函数积分\n"));
	builtin_functions.insert(pair<string, function_type>("integrate", function_type(integrate, 3)));

	command_instruction.insert(map < string, string >::value_type("diff", "diff(函数名,求导数值的位置)：求函数导数\n"));
	builtin_functions.insert(pair<string, function_type>("diff", function_type(diff, 2)));

	command_instruction.insert(map < string, string >::value_type("cos", "cos(x)：\n"));
	builtin_functions.insert(pair<string, function_type>("cos", function_type(cos, 1)));

	command_instruction.insert(map < string, string >::value_type("sin", "sin(x)：\n"));
	builtin_functions.insert(pair<string, function_type>("sin", function_type(sin, 1)));

	command_instruction.insert(map < string, string >::value_type("tan", "tan(x)：\n"));
	builtin_functions.insert(pair<string, function_type>("tan", function_type(tan, 1)));

	command_instruction.insert(map < string, string >::value_type("acos", "acos(x)：\n"));
	builtin_functions.insert(pair<string, function_type>("acos", function_type(acos, 1)));

	command_instruction.insert(map < string, string >::value_type("asin", "asin(x)：\n"));
	builtin_functions.insert(pair<string, function_type>("asin", function_type(asin, 1)));

	command_instruction.insert(map < string, string >::value_type("atan", "atan(x)：\n"));
	builtin_functions.insert(pair<string, function_type>("atan", function_type(atan, 1)));
}

void operator_pow(std::vector<std::string> *params) {
	using namespace std;
	int n = params->size();
	if (n < 2) {
		throw  string("表达式错误");
	}
	if (symbols.find((*params)[n - 1]) != symbols.end()) {
		(*params)[n - 1] = to_string(*symbols[(*params)[n - 1]]);
	}
	if (symbols.find((*params)[n - 2]) != symbols.end()) {
		(*params)[n - 2] = to_string(*symbols[(*params)[n - 2]]);
	}
	(*params)[n - 2] = to_string(pow(stod((*params)[n - 2]), stod((*params)[n - 1])));
	params->pop_back();
}
void operator_mutiply(std::vector<std::string> *params) {
	using namespace std;
	int n = params->size();
	if (n < 2) {
		throw  string("表达式错误");
	}
	if (result_matrix.size()==0 && matrixes.find((*params)[n - 1]) != matrixes.end() && matrixes.find((*params)[n - 2]) != matrixes.end()) {
		result_matrix.push_back(*matrixes[(*params)[n - 2]] * *matrixes[(*params)[n - 1]]);
		params->pop_back();
		params->back() = "(";
	}
	else if (result_matrix.size() != 0) {
		shared_ptr<Matrix<double>> left, right;
		int m = result_matrix.size();
		if ((*params)[n - 2] == "(" && (*params)[n - 1] == "(") {
			left = result_matrix[m - 2];
			right = result_matrix[m - 1];
			result_matrix.pop_back();
			result_matrix.pop_back();
		}
		else if ((*params)[n - 2] != "(" && (*params)[n - 1] != "(") {
			left = matrixes[(*params)[n - 2]];
			right = matrixes[(*params)[n - 1]];
		}
		else if ((*params)[n - 2] == "(" && (*params)[n - 1] != "(") {
			left = result_matrix[m - 1];
			right = matrixes[(*params)[n - 1]];
			result_matrix.pop_back();
		}
		else {
			left = matrixes[(*params)[n - 2]];
			right = result_matrix[m - 1];
			result_matrix.pop_back();
		}
		params->pop_back();
		params->back() = "(";
		result_matrix.push_back(*left * *right);
	}
	else {
		if (symbols.find((*params)[n - 1]) != symbols.end()) {
			(*params)[n - 1] = to_string(*symbols[(*params)[n - 1]]);
		}
		if (symbols.find((*params)[n - 2]) != symbols.end()) {
			(*params)[n - 2] = to_string(*symbols[(*params)[n - 2]]);
		}
		(*params)[n - 2] = to_string(stod((*params)[n - 2])*stod((*params)[n - 1]));
		params->pop_back();
	}

}
void operator_divide(std::vector<std::string> *params) {
	using namespace std;
	int n = params->size();
	if (n < 2) {
		throw  string("表达式错误");
	}
	if (stod((*params)[n - 1])==0) {
		throw  string("表达式除以0");
	}
	if (symbols.find((*params)[n - 1]) != symbols.end()) {
		(*params)[n - 1] = to_string(*symbols[(*params)[n - 1]]);
	}
	if (symbols.find((*params)[n - 2]) != symbols.end()) {
		(*params)[n - 2] = to_string(*symbols[(*params)[n - 2]]);
	}
	(*params)[n - 2] = to_string(stod((*params)[n - 2]) / stod((*params)[n - 1]));
	params->pop_back();
}
void operator_add(std::vector<std::string> *params) {
	using namespace std;
	int n = params->size();
	if (n < 2) {
		throw  string("表达式错误");
	}
	if (result_matrix.size() == 0 && matrixes.find((*params)[n - 1]) != matrixes.end() && matrixes.find((*params)[n - 2]) != matrixes.end()) {
		result_matrix.push_back(*matrixes[(*params)[n - 2]] + *matrixes[(*params)[n - 1]]);
		params->pop_back();
		params->back() = "(";
	}
	else if (result_matrix.size() != 0) {
		shared_ptr<Matrix<double>> left, right;
		int m = result_matrix.size();
		if ((*params)[n - 2] == "(" && (*params)[n - 1] == "(") {
			left = result_matrix[m - 2];
			right = result_matrix[m - 1];
			result_matrix.pop_back();
			result_matrix.pop_back();
		}
		else if ((*params)[n - 2] != "(" && (*params)[n - 1] != "(") {
			left = matrixes[(*params)[n - 2]];
			right = matrixes[(*params)[n - 1]];
		}
		else if ((*params)[n - 2] == "(" && (*params)[n - 1] != "(") {
			left = result_matrix[m - 1];
			right = matrixes[(*params)[n - 1]];
			result_matrix.pop_back();
		}
		else {
			left = matrixes[(*params)[n - 2]];
			right = result_matrix[m - 1];
			result_matrix.pop_back();
		}
		params->pop_back();
		params->back() = "(";
		result_matrix.push_back(*left + *right);
	}
	else {
		if (symbols.find((*params)[n - 1]) != symbols.end()) {
			(*params)[n - 1] = to_string(*symbols[(*params)[n - 1]]);
		}
		if (symbols.find((*params)[n - 2]) != symbols.end()) {
			(*params)[n - 2] = to_string(*symbols[(*params)[n - 2]]);
		}
		(*params)[n - 2] = to_string(stod((*params)[n - 2]) + stod((*params)[n - 1]));
		params->pop_back();
	}
}
void operator_minus(std::vector<std::string> *params) {
	using namespace std;
	int n = params->size();
	if (n < 2) {
		throw  string("表达式错误");
	}
	if (result_matrix.size() == 0 && matrixes.find((*params)[n - 1]) != matrixes.end() && matrixes.find((*params)[n - 2]) != matrixes.end()) {
		result_matrix.push_back(*matrixes[(*params)[n - 2]] - *matrixes[(*params)[n - 1]]);
		params->pop_back();
		params->back() = "(";
	}
	else if (result_matrix.size() != 0) {
		shared_ptr<Matrix<double>> left, right;
		int m = result_matrix.size();
		if ((*params)[n - 2] == "(" && (*params)[n - 1] == "(") {
			left = result_matrix[m - 2];
			right = result_matrix[m - 1];
			result_matrix.pop_back();
			result_matrix.pop_back();
		}
		else if ((*params)[n - 2] != "(" && (*params)[n - 1] != "(") {
			left = matrixes[(*params)[n - 2]];
			right = matrixes[(*params)[n - 1]];
		}
		else if ((*params)[n - 2] == "(" && (*params)[n - 1] != "(") {
			left = result_matrix[m - 1];
			right = matrixes[(*params)[n - 1]];
			result_matrix.pop_back();
		}
		else {
			left = matrixes[(*params)[n - 2]];
			right = result_matrix[m - 1];
			result_matrix.pop_back();
		}
		params->pop_back();
		params->back() = "(";
		result_matrix.push_back(*left - *right);
	}
	else {
		if (symbols.find((*params)[n - 1]) != symbols.end()) {
			(*params)[n - 1] = to_string(*symbols[(*params)[n - 1]]);
		}
		if (symbols.find((*params)[n - 2]) != symbols.end()) {
			(*params)[n - 2] = to_string(*symbols[(*params)[n - 2]]);
		}
		(*params)[n - 2] = to_string(stod((*params)[n - 2]) - stod((*params)[n - 1]));
		params->pop_back();
	}

}
void operator_assign() {
	using namespace std;
	if (var_name == "") return;
	if (result_matrix.size() == 1) {
		if (matrixes.find(var_name) == matrixes.end()) {
			shared_ptr<Matrix<double>> new_matrix(new Matrix<double>());
			matrixes.insert(map < string, shared_ptr<Matrix<double>> >::value_type(var_name, new_matrix));
		}
		*(matrixes[var_name]) = *result_matrix[0];
	}
	else if (function_params.size() == 1 && matrixes.find(function_params.back()) != matrixes.end()) {
		if (matrixes.find(var_name) == matrixes.end()) {
			shared_ptr<Matrix<double>> new_matrix(new Matrix<double>());
			matrixes.insert(map < string, shared_ptr<Matrix<double>> >::value_type(var_name, new_matrix));
		}
		*(matrixes[var_name]) = *(matrixes[function_params.back()]);
	}
	else if (function_params.size() == 1) {
		if (symbols.find(function_params.back()) != symbols.end()) {
			function_params[0] = to_string(*symbols[function_params[0]]);
		}
		if (symbols.find(var_name) == symbols.end()) {
			symbols.insert(map < string, shared_ptr<double>>::value_type(var_name, nullptr));
		}
		if (symbols[var_name] == nullptr) {
			symbols[var_name] = shared_ptr<double>(new double);
		}
		*symbols[var_name] = stod(function_params.back());
	}

}
void save_name(const std::string &name) {
	var_name = name;
}

boost::spirit::qi::rule<std::string::iterator, std::string()>value;
boost::spirit::qi::rule<std::string::iterator, std::string()> variable_name;
boost::spirit::qi::rule<std::string::iterator, std::string()> function_name_rule;
boost::spirit::qi::rule<std::string::iterator, std::string()> command_name_rule;
boost::spirit::qi::rule<std::string::iterator, std::string()>command_rule;//匹配基本指令
boost::spirit::qi::rule<std::string::iterator, boost::spirit::standard::space_type, std::string()> function_rule;//匹配自定义函数
boost::spirit::qi::rule<std::string::iterator, boost::spirit::standard::space_type, std::string()> priority0,//函数运算、读取数值和括号操作
priority1,//幂运算
priority2,//乘除法运算
priority3,//加减法运算
assign_rule,
priority4;//赋值

void init_parsing() {
	using namespace boost::spirit::qi;
	using boost::phoenix::push_back;
	using boost::bind;

	variable_name = +(char_("a", "z")[_val += _1] | char_("A", "Z")[_val += _1] | char_("0", "9")[_val += _1] | char_("_")[_val += _1]);
	command_name_rule = +(char_("a", "z")[_val += _1] | char_("A", "Z")[_val += _1] | char_("0", "9")[_val += _1]);
	command_rule = command_name_rule[boost::phoenix::ref(command_name) = _1] >> *(' ' >> variable_name[push_back(boost::phoenix::ref(command_params), _1)]);
	
	value = -(char_('-')[_val += _1]) >> +(char_('.')[_val += _1] | char_("0", "9")[_val += _1]);
	function_name_rule = +(char_("a", "z")[_val += _1] | char_("A", "Z")[_val += _1] | char_("0", "9")[_val += _1]) >> *(char_(' ')) >> char_('(')[push_back(boost::phoenix::ref(function_params), "(")];
	function_rule = function_name_rule[push_back(boost::phoenix::ref(function_names), _1)] >> -(priority3 >> *(',' >> priority3)) >> char_(')')[bind(run_function, &function_names, &function_params)];

	priority0 = function_rule | variable_name[push_back(boost::phoenix::ref(function_params), _1)] | value[push_back(boost::phoenix::ref(function_params), _1)] |  ('(' >> priority3 >> ')');
	priority1 = priority0 >> *('^' >> priority0[bind(operator_pow, &function_params)]);
	priority2 = priority1 >> *(('*' >> priority1[bind(operator_mutiply, &function_params)]) | ('/' >> priority1[bind(operator_divide, &function_params)]));
	priority3 = priority2 >> *(('+' >> priority2[bind(operator_add, &function_params)]) | ('-' >> priority2[bind(operator_minus, &function_params)]));
	assign_rule = variable_name[_val = _1] >> '=';
	priority4 = -(assign_rule[bind(save_name, boost::lambda::_1)]) >> priority3[bind(operator_assign)];
}

bool command_parsing(std::string &instruction) {
	using namespace boost::spirit::qi;
	using boost::phoenix::ref;
	using boost::phoenix::push_back;
	using boost::bind;
	auto begin = instruction.begin();
	auto end = instruction.end();
	bool succeed = false;
	parse(begin, end, command_rule);//指令匹配
	if (begin== end) {
		succeed=run_command(command_name, &command_params);
	}
	command_params.clear();
	if (!succeed) {//函数匹配
		begin = instruction.begin();
		try {
			phrase_parse(begin, end, priority4, space);
			succeed = true;
			if (begin == end) {
				if (function_params.size() == 1) {
					if (result_matrix.size() == 1)result_matrix[0]->print();
					else std::cout << function_params[0] << std::endl;
				}
				else {
					succeed = false;
					std::cout << "参数错误" << std::endl;
				}
			}
			else {
				succeed = false;
				std::cout << "表达式错误" << std::endl;
			}
		}
		catch (const std::string error) {
			std::cout << error << std::endl;
		}
	}
	function_names.clear();
	function_params.clear();
	result_matrix.clear();
	var_name = "";
	return succeed;
}