#include "statement.h"

namespace client {
	namespace parser {

		statement::statement(): statement::base_type(statement_list), expr() {

			qi::raw_type raw;
			qi::lexeme_type lexeme;
			qi::alpha_type alpha;
			qi::alnum_type alnum;
			qi::double_type double_;

			statement_list =
				+(
					 variable_declaration
					| assignment
					| (expr >> ';')
					)
				;


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