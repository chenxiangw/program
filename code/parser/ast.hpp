#pragma once

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/variant/recursive_variant.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/optional.hpp>
#include <list>

namespace client { 
	///////////////////////////////////////////////////////////////////////////
	//  The AST
	///////////////////////////////////////////////////////////////////////////
	namespace ast{

		struct nil {};
		struct signed_;
		struct function_call;
		struct expression;

		struct variable {
			std::string var;
		};

		struct vector {
			std::list<double> values;
		};

		struct matrix {
			std::list<vector> rhs;
		};

		typedef boost::variant<
			nil
			, double
			, variable
			, matrix
			, vector
			, boost::recursive_wrapper<signed_>
			, boost::recursive_wrapper<function_call>
			, boost::recursive_wrapper<expression>
		>operand;

		struct signed_{
			char sign;
			operand operand_;
		};

		struct operation{
			char operator_;
			operand operand_;
		};

		typedef boost::variant<
			double
			, variable
		>parameter_type;

		struct parameter {
			parameter_type param;
		};

		struct function_call{
			std::string function_name;
			std::list<parameter> args;
		};

		struct expression{
			operand first;
			std::list<operation> rest;
		};

		struct assignment{
			std::string lhs;
			expression rhs;
		};

		struct variable_declaration{
			std::list<std::string> vars;
			std::list<expression> rhs;
		};

		typedef boost::variant<
			variable_declaration
			, assignment
			, expression
			>
		statement;

		typedef std::list<statement> statement_list;

	}
}

BOOST_FUSION_ADAPT_STRUCT(
	client::ast::variable,
	(std::string, var)
)

BOOST_FUSION_ADAPT_STRUCT(
	client::ast::vector,
	(std::list<double>, values)
)

BOOST_FUSION_ADAPT_STRUCT(
	client::ast::matrix,
	(std::list<client::ast::vector>, rhs)
)

BOOST_FUSION_ADAPT_STRUCT(
    client::ast::signed_,
    (char, sign)
    (client::ast::operand, operand_)
)

BOOST_FUSION_ADAPT_STRUCT(
    client::ast::operation,
    (char, operator_)
    (client::ast::operand, operand_)
)

BOOST_FUSION_ADAPT_STRUCT(
	client::ast::parameter,
	(client::ast::parameter_type, param)
)

BOOST_FUSION_ADAPT_STRUCT(
	client::ast::function_call,
	(std::string, function_name)
	(std::list<client::ast::parameter>, args)
)


BOOST_FUSION_ADAPT_STRUCT(
    client::ast::expression,
    (client::ast::operand, first)
    (std::list<client::ast::operation>, rest)
)

BOOST_FUSION_ADAPT_STRUCT(
	client::ast::assignment,
	(std::string, lhs)
	(client::ast::expression, rhs)
)

BOOST_FUSION_ADAPT_STRUCT(
	client::ast::variable_declaration,
	(std::list<std::string>, vars)
	(std::list<client::ast::expression>, rhs)
)
