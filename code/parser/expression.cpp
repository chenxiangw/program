

#include "expression.h"
#include <boost/spirit/include/phoenix_function.hpp>

namespace client {
	namespace parser {
		typedef std::string::const_iterator iterator_type;
		expression::expression(): expression::base_type(expr) {

			qi::char_type char_;
			qi::raw_type raw;
			qi::lexeme_type lexeme;
			qi::alpha_type alpha;
			qi::alnum_type alnum;
			qi::double_type double_;

			expr =
				power_expr.alias()
				;

			power_expr =
				additive_expr
				>> *(char_('^') > additive_expr)
				;

			additive_expr =
				multiplicative_expr
				>> *((char_('+') > multiplicative_expr)
					| (char_('-') > multiplicative_expr)
					)
				;

			multiplicative_expr =
				unary_expr
				>> *((char_('*') > unary_expr)
					| (char_('/') > unary_expr)
					)
				;

			unary_expr =
				primary_expr
				| (char_('-') > primary_expr)
				| (char_('+') > primary_expr)
				;

			primary_expr =
				double_
				| matrix
				| vector
				| function_call
				| variable
				| '(' > expr > ')'
				;

			matrix =
				("[" >> (vector%',') >> ']');
				

			vector = '[' >> (double_%',') >> ']';

			function_call =
				(function_name >> '(')
				>> -(parameter % ',')
				>>   ')'
				;

			parameter = (double_ | variable);

			variable = function_name;

			function_name =
				raw[lexeme[(alpha | '_') >> *(alnum | '_')]]
				;

		}
	}
}