#pragma once


#include <boost/spirit/include/qi.hpp>
#include "ast.hpp"
#include <vector>

namespace client { 
	namespace parser{
		namespace qi = boost::spirit::qi;
		namespace ascii = boost::spirit::ascii;

		///////////////////////////////////////////////////////////////////////////////
		//  The expression grammar
		///////////////////////////////////////////////////////////////////////////////

		typedef std::string::const_iterator iterator_type;
		struct expression : qi::grammar<iterator_type, ast::expression(), ascii::space_type>{
			expression();
			qi::rule<iterator_type, ast::expression(), ascii::space_type>
				expr, additive_expr, multiplicative_expr, power_expr
				;

			qi::rule<iterator_type, ast::function_call(), ascii::space_type >function_call;

			qi::rule<iterator_type, ast::matrix(), ascii::space_type> matrix;
			qi::rule<iterator_type, ast::vector(), ascii::space_type> vector;

			qi::rule<iterator_type, ast::operand(), ascii::space_type>
				unary_expr, primary_expr
				;

			qi::rule<iterator_type, std::string(), ascii::space_type> function_name;
			qi::rule<iterator_type, ast::variable(), ascii::space_type> variable;
			qi::rule<iterator_type, ast::parameter_type(), ascii::space_type> parameter;
		};
	}
}


