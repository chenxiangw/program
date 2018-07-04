#pragma once

#include "ast.hpp"
#include "program.hpp"
#include "statement.h"
#include"../Global_variable.hpp"
#include <vector>
#include <map>

namespace client { 


	///////////////////////////////////////////////////////////////////////////
	//  The Compiler
	///////////////////////////////////////////////////////////////////////////
	struct compiler {
		typedef bool result_type;

		compiler(client::program& program, Global_variable &global_variable)
			: program(program), global_variable(global_variable) {
		}

		bool operator()(ast::nil) const { BOOST_ASSERT(0); return false; }
		bool operator()(double x) const;
		bool operator()(ast::variable const& x) const;
		bool operator()(ast::vector const& x) const;
		bool operator()(ast::matrix const& x) const;
		bool operator()(ast::operation const& x) const;
		bool operator()(ast::signed_ const& x) const;
		bool operator()(ast::function_call const& x)const;
		bool operator()(ast::expression const& x) const;
		bool operator()(ast::assignment const& x) const;
		bool operator()(ast::variable_declaration const& x) const;
		bool operator()(ast::define_function const& x) const;
		bool operator()(ast::statement_list const& x) const;

		client::program& program;
		Global_variable &global_variable;

	};

	///////////////////////////////////////////////////////////////////////////
	//  The user_function_compiler
	///////////////////////////////////////////////////////////////////////////
	struct user_function_compiler {
		typedef bool result_type;

		user_function_compiler(client::function_program & program)
			: program(program) {
		}

		bool operator()(ast::nil) const { BOOST_ASSERT(0); return false; }
		bool operator()(double x) const;
		bool operator()(ast::variable const& x) const;
		bool operator()(ast::vector const& x) const;
		bool operator()(ast::matrix const& x) const;
		bool operator()(ast::operation const& x) const;
		bool operator()(ast::signed_ const& x) const;
		bool operator()(ast::function_call const& x)const;
		bool operator()(ast::expression const& x) const;
		bool operator()(ast::assignment const& x) const;
		bool operator()(ast::variable_declaration const& x) const;
		bool operator()(ast::define_function const& x) const;
		bool operator()(ast::statement_list const& x) const;

		client::function_program & program;
	};
}

