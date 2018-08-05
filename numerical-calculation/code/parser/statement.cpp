#include "statement.h"

namespace client {
	namespace parser {

		statement::statement(): statement::base_type(statement_list), expr() {

			qi::raw_type raw;
			qi::lexeme_type lexeme;
			qi::alpha_type alpha;
			qi::alnum_type alnum;
			qi::double_type double_;
			qi::char_type char_;

			statement_list =
				+(
					define_function
					| variable_declaration
					| assignment
					| (expr >> ';')
					)
				;

			define_function =
				lexeme["function"]
				>> function_info
				>> content
				;

			content = '{' >> *(char_ - "};") >> "};";
			
			function_info=
				'[' >> -(output%',') >> "]="
				>> function_name
				>> '(' >> -(input%',') >> ')'
				;

			function_name = raw[lexeme[(alpha | '_') >> *(alnum | '_')]];
			output = param;
			input = type>>' '>>param;
			type = raw[lexeme[(alpha | '_') >> *(alnum | '_')]];
			param = raw[lexeme[(alpha | '_') >> *(alnum | '_')]];

			variable_declaration =
				lexeme["var" >> !(alnum | '_')]
				>> ( declare_variable % ',')
				>> -('=' >> expr)
				>> ';'
				;

			declare_variable =
				raw[lexeme[(alpha | '_') >> *(alnum | '_')]]
				;

			assignment =
				assign_variable
				>> '='
				>> expr
				>> ';'
				;

			assign_variable =
				raw[lexeme[(alpha | '_') >> *(alnum | '_')]]
				;
		}
	}
}