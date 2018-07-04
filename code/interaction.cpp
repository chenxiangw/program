#include<iostream>
#include<vector>
#include<string>
#include<set>

#include "parser/vm.hpp"


void interaction() {
	Global_variable global_variable;
	global_variable.initiate();
	client::program program;
	client::compiler compile(program, global_variable);
	client::vmachine vm(global_variable);

	std::string source;
	while (true) {
		source = "";
		std::cout << "ÇëÊäÈëÖ¸Áî£º" << std::endl;
		std::string line = "";
		while (std::getline(std::cin, line)){
			source += line;
			if (line.size() > 0 && line.back() == ';')break;
		}

		typedef std::string::const_iterator iterator_type;
		iterator_type iter = source.begin();
		iterator_type end = source.end();
		client::parser::statement parser;
		boost::spirit::ascii::space_type space;
		client::ast::statement_list ast;
		bool success = phrase_parse(iter, end, parser, space, ast);//½âÎöÓï·¨Ê÷

		std::cout << "-------------------------\n";

		if (success && iter == end) {
			if (compile(ast)) {//±àÒë
				program.print_assembler();
				std::cout << "compile success\n";
				std::cout << "-------------------------\n";
				if (vm.execute(program))vm.print_stack();
				else {
					std::cout << "execute failure\n";
				}
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