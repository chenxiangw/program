#pragma once

#include "expression.h"

namespace client {
	namespace parser{
		///////////////////////////////////////////////////////////////////////////////
		//  The statement grammar
		///////////////////////////////////////////////////////////////////////////////
		typedef std::string::const_iterator iterator_type;
		struct statement : qi::grammar<iterator_type, ast::statement_list(), ascii::space_type>{
			statement();
			expression expr;
			qi::rule<iterator_type, ast::statement_list(), ascii::space_type> statement_list;
			qi::rule<iterator_type, ast::variable_declaration(), ascii::space_type> variable_declaration;
			qi::rule<iterator_type, ast::assignment(), ascii::space_type> assignment;
			qi::rule<iterator_type, std::string(), ascii::space_type>
				assign_variable, declare_variable;
		};
	}
}



