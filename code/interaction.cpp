#include<iostream>
#include<vector>
#include<string>

#include "./parser/statement.h"
#include "./parser/vm.hpp"
#include "./parser/compiler.h"

#include"builtin_functions.h"
#include"./core/Matrix.hpp"
#include "./core/Matrix.hpp"

namespace user_input {
	std::map<std::string, std::pair<bool, double>> variables;
	std::map < std::string, std::shared_ptr<Matrix<double>>> matrixes;
}

class Add_builtin_function {
public:
	Add_builtin_function(std::map < std::string, void(*)(client::vmachine  *vm) > &builtin_functions, std::map < std::string, std::vector<std::string>  > &function_exsit)
		:builtin_functions(builtin_functions), function_exsit(function_exsit) {
	}
	void add(std::string name, std::vector<std::string> param_type, void(*address)(client::vmachine  *vm)) {
		builtin_functions.insert(std::pair <std::string, void(*)(client::vmachine *vm)>(name, address));
		function_exsit.insert(std::pair<std::string, std::vector<std::string> >(name, param_type));
	}
private:
	std::map < std::string, void(*)(client::vmachine  *vm) > &builtin_functions;
	std::map < std::string, std::vector<std::string>  > &function_exsit;
};

void init_program(std::map < std::string, void(*)(client::vmachine  *vm) > &builtin_functions, std::map < std::string, std::vector<std::string>  > &function_exsit) {
	class Add_builtin_function add_builtin_function(builtin_functions, function_exsit);
	std::string double_type = typeid(double).name();
	std::string matrix_type = typeid(std::shared_ptr<Matrix<double>>).name();
	std::string string_type = typeid(std::string).name();

	add_builtin_function.add("cos", { double_type }, cos_builtin);
	add_builtin_function.add("sin", { double_type }, sin_builtin);
	add_builtin_function.add("tan", { double_type }, tan_builtin);
	add_builtin_function.add("acos", { double_type }, acos_builtin);
	add_builtin_function.add("asin", { double_type }, asin_builtin);
	add_builtin_function.add("atan", { double_type }, atan_builtin);

	add_builtin_function.add("gamma", { double_type }, gamma_builtin);
	add_builtin_function.add("ln_gamma", { double_type }, ln_gamma_builtin);
	add_builtin_function.add("gamma_P", { double_type, double_type }, gamma_P_builtin);
	add_builtin_function.add("gamma_Q", { double_type, double_type }, gamma_Q_builtin);
	add_builtin_function.add("factorial", { double_type }, factorial_builtin);
	add_builtin_function.add("binomial", { double_type, double_type }, binomial_builtin);
	add_builtin_function.add("beta", { double_type, double_type }, beta_builtin);
	add_builtin_function.add("incomplete_beta", { double_type, double_type, double_type }, incomplete_beta_builtin);
	add_builtin_function.add("error_function", { double_type }, error_function_builtin);
	add_builtin_function.add("Legendre_p", { double_type, double_type }, Legendre_p_builtin);
	add_builtin_function.add("Legendre_q", { double_type, double_type }, Legendre_q_builtin);

	add_builtin_function.add("inverse", { matrix_type }, inverse_builtin);
	add_builtin_function.add("determinant", { matrix_type }, determinant_builtin);
	add_builtin_function.add("solve", { matrix_type, matrix_type }, solve_builtin);
}

void interaction() {
	std::map < std::string, void(*)(client::vmachine  *vm) > builtin_functions;
	std::map < std::string, std::vector<std::string> > function_exsit;
	init_program(builtin_functions, function_exsit);

	client::code_gen::program program(user_input::variables, user_input::matrixes , function_exsit);
	client::code_gen::compiler compile(program);
	client::vmachine vm(user_input::variables, user_input::matrixes, builtin_functions);

	std::string source;
	while (true) {
		source = "";
		std::cout << "请输入指令（quit退出）：" << std::endl;
		std::getline(std::cin, source);
		if (source == "quit") break;

		typedef std::string::const_iterator iterator_type;
		iterator_type iter = source.begin();
		iterator_type end = source.end();
		client::parser::statement parser;
		boost::spirit::ascii::space_type space;
		client::ast::statement_list ast;
		bool success = phrase_parse(iter, end, parser, space, ast);//解析语法树

		std::cout << "-------------------------\n";

		if (success && iter == end) {
			if (compile(ast)) {//编译
				program.print_assembler();
				std::cout << "Success\n";
				std::cout << "-------------------------\n";
				if (vm.execute(program))
					vm.print_stack();//运行
			}
			else {
				std::cout << "Compile failure\n";
			}
		}
		else {
			std::cout << "Parse failure\n";
		}

		program.clear();
		vm.clear();
	}

}