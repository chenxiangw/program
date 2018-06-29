
#include "./parser/statement.hpp"
#include "./parser/vm.hpp"
#include "./parser/compiler.hpp"


namespace saved_input {
	std::map<std::string, std::pair<bool, double>> variables;
	typedef  std::pair<void(*)(const std::vector<std::string> *), int> function_type;
	std::map < std::string, function_type > builtin_functions;
}

///////////////////////////////////////////////////////////////////////////////
//  Main program
///////////////////////////////////////////////////////////////////////////////
int main(){

	client::vmachine vm;                      
	client::code_gen::program program(saved_input::variables, saved_input::builtin_functions);

    std::string source;
	while (std::getline(std::cin, source)) {
		if (source.empty())
			break;

		typedef std::string::const_iterator iterator_type;
		iterator_type iter = source.begin();
		iterator_type end = source.end();

		client::error_handler<iterator_type> error_handler(iter, end);     
		client::parser::statement parser(error_handler);
		client::code_gen::compiler compile(program, error_handler);       
		boost::spirit::ascii::space_type space;
		client::ast::statement_list ast;
		bool success = phrase_parse(iter, end, parser, space, ast);

		std::cout << "-------------------------\n";

		if (success && iter == end){
			if (compile(ast)){
				program.print_assembler();
				std::cout << "Success\n";
				std::cout << "-------------------------\n";
				vm.execute(program);

			}
			else{
				std::cout << "Compile failure\n";
			}
		}
		else{
			std::cout << "Parse failure\n";
		}

		std::cout << "-------------------------\n\n";
		program.clear();
	}
    return 0;
}


