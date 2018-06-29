#pragma once

#include "ast.hpp"
#include "../data_type.h"
#include <vector>
#include <map>
#include <boost/function.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_function.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

namespace client { 
	namespace code_gen{
		///////////////////////////////////////////////////////////////////////////
		//  The Program
		///////////////////////////////////////////////////////////////////////////
		struct program{
			program(std::map<std::string, std::pair<bool, double>> &variables,
				std::map < std::string, std::shared_ptr<Matrix<double>>> &matrixes,
				std::map < std::string, std::vector<std::string> > &function_calls)
				:variables(variables), function_calls(function_calls), matrixes(matrixes){
			}

			void save_code(int code);
			void save_value(double value);
			void save_variable_name(std::string const &name);
			void save_matrix(std::shared_ptr<Matrix<double>> &matrix);

			int& operator[](std::size_t i) { return code[i]; }
			int const& operator[](std::size_t i) const { return code[i]; }

			std::pair<bool, double> const* find_var(std::string const& name) const;
			void add_var(std::string const& name);
			void set_var(std::string const& name,double value);

			std::vector<std::string> const* find_func(std::string const& name) const;
			std::shared_ptr<Matrix<double>> const * find_matrix(std::string const& name)const;
			void add_matrix(std::string const& name, std::shared_ptr<Matrix<double>> &new_matrix);

			std::vector<int> const& get_code() const { return code; }
			std::vector<double> const& get_value() const { return value; }
			std::vector<std::string> const& get_variables_loc() const { return variable_names; }
			std::vector< std::shared_ptr<Matrix<double>>> const& get_matrix_loc() const { return matrix_loc; }

			void clear() {
				code.clear();
				value.clear();
				matrix_loc.clear();
				variable_names.clear();
			}
			void print_assembler() const;

		private:

			
			std::map < std::string, std::vector<std::string> > &function_calls;
			std::map < std::string, std::shared_ptr<Matrix<double>>> &matrixes;
			std::map<std::string, std::pair<bool, double>> &variables;
			std::vector<int> code;
			std::vector<double> value;
			std::vector< std::shared_ptr<Matrix<double>>> matrix_loc;
			std::vector<std::string> variable_names;
		};

		///////////////////////////////////////////////////////////////////////////
		//  The Compiler
		///////////////////////////////////////////////////////////////////////////
		struct compiler{
			typedef bool result_type;


			compiler(client::code_gen::program& program)
			  : program(program){
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
			bool operator()(ast::statement_list const& x) const;

			client::code_gen::program& program;


		};
	}
}

