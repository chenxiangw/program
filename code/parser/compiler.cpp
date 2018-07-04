
#include "compiler.h"
#include "vm.hpp"
#include <boost/foreach.hpp>
#include <boost/variant/apply_visitor.hpp>


namespace client { 
	std::map<std::string, std::string> type_name = { {"double",typeid(std::shared_ptr<double>).name()},
	{"matrix",typeid(std::shared_ptr<Matrix<double>>).name()},
	{"function",typeid(std::shared_ptr<function_program>).name()},
	};
	///////////////////////////////////////////////////////////////////////////
	//  The Compiler
	///////////////////////////////////////////////////////////////////////////
	bool compiler::operator()(double x) const{
		program.save_code(op_data);
		program.save_data(x);
		return true;
	}

	bool compiler::operator()(ast::variable const& x) const {
		program.save_code(op_load);
		program.save_variable_name(x.var);
		return true;
	}

	bool compiler::operator()(ast::vector const& x) const {
		size_t n = x.values.size();
		std::shared_ptr<Matrix<double>> new_matrix;
		if (n != 0) {
			std::list<double>::const_iterator val = x.values.begin();
			new_matrix = std::shared_ptr<Matrix<double>>(new Matrix<double>(1, n));
			for (int column = 0; val != x.values.end(); ++val, ++column) {
				(*new_matrix)[0][column] = *val;
			}
		}
		else {
			new_matrix = std::shared_ptr<Matrix<double>>(new Matrix<double>());
		}
		program.save_code(op_data);
		program.save_data(new_matrix);
		return true;
	}

	bool compiler::operator()(ast::matrix const& x) const {
		size_t m = x.rhs.size();
		std::shared_ptr<Matrix<double>> new_matrix;
		if (m != 0) {
			std::list<client::ast::vector>::const_iterator vec = x.rhs.begin();
			size_t n = vec->values.size();
			new_matrix = std::shared_ptr<Matrix<double>>(new Matrix<double>(m, n));
			for (int row = 0; vec != x.rhs.end(); ++vec, ++row) {
				if (vec->values.size() != n) {
					std::cout << "矩阵列数不一致" << std::endl;
					return false;
				}
				std::list<double>::const_iterator val = vec->values.begin();
				for (int column = 0; val != vec->values.end(); ++val, ++column) {
					(*new_matrix)[row][column] = *val;
				}
			}
		}
		else {
			new_matrix = std::shared_ptr<Matrix<double>>(new Matrix<double>());
		}
		program.save_code(op_data);
		program.save_data(new_matrix);
		return true;
	}

	bool compiler::operator()(ast::operation const& x) const{
		if (!boost::apply_visitor(*this, x.operand_))
			return false;
		switch (x.operator_){
			case '+': program.save_code(op_add); break;
			case '-': program.save_code(op_sub); break;
			case '*': program.save_code(op_mul); break;
			case '/': program.save_code(op_div); break;
			case '^': program.save_code(op_pow); break;
			default: BOOST_ASSERT(0); return false;
		}
		return true;
	}

	bool compiler::operator()(ast::signed_ const& x) const{
		if (!boost::apply_visitor(*this, x.operand_))
			return false;
		switch (x.sign){
			case '-': program.save_code(op_neg); break;
			case '+': break;
			default: BOOST_ASSERT(0); return false;
		}
		return true;
	}

	bool compiler::operator()(ast::function_call const& x)const {
		program.save_code(op_func_call);
		program.save_variable_name(x.function_name);
		BOOST_FOREACH(ast::parameter_type const& param, x.args) {
			program.save_code(op_param);
			if (!boost::apply_visitor(*this, param)) {
				return false;
			}
		}
		return true;
	}

	bool compiler::operator()(ast::expression const& x) const{
		if (!boost::apply_visitor(*this, x.first))
			return false;
		BOOST_FOREACH(ast::operation const& oper, x.rest){
			if (!(*this)(oper))
				return false;
		}
		return true;
	}

	bool compiler::operator()(ast::assignment const& x) const{
		if (!(*this)(x.rhs))return false;
		program.save_code(op_store);
		program.save_variable_name(x.lhs);
		return true;
	}

	bool compiler::operator()(ast::variable_declaration const& x) const {
		if (x.rhs.size()>0) {//存在赋值语句
			BOOST_FOREACH(ast::expression const& expr, x.rhs){
				if (!(*this)(expr))
					return false;
			}
			for (auto &elem : x.vars) {
				program.save_code(op_store_var);
				program.save_variable_name(elem);
			}
		}
		else {
			for (auto &elem : x.vars) {
				global_variable.declara_var(elem);
			}
		}
		return true;
	}

	bool compiler::operator()(ast::define_function const& x) const {
		std::vector<std::pair<std::string, std::string >> input;
		std::vector<std::string> output;
		for (auto &elem : x.info.in) {
			auto p = type_name.find(elem.type);
			if (p == type_name.end())return false;
			input.push_back(std::pair<std::string, std::string >(elem.param, p->second));
		}
		for (auto &elem : x.info.out) {
			output.push_back(elem.param);
		}
		std::shared_ptr<client::function_program> new_function(new function_program);
		new_function->initiate(x.content, output, input, global_variable.get_illegal_names_sys());
		typedef std::string::const_iterator iterator_type;
		iterator_type iter = x.content.begin();
		iterator_type end = x.content.end();
		client::parser::statement parser;
		boost::spirit::ascii::space_type space;
		client::ast::statement_list ast;
		bool success = phrase_parse(iter, end, parser, space, ast);//解析语法树
		if (!success || iter != end)return false;
		client::user_function_compiler compiler(*new_function);
		if (!compiler(ast))return false;

		global_variable.assign_function(x.info.function_name, new_function);


		return true;
	}

	bool compiler::operator()(ast::statement_list const& x) const{
		BOOST_FOREACH(ast::statement const& statement, x){
			if (!boost::apply_visitor(*this, statement)){
				return false;
			}
		}
		return true;
	}
	
	///////////////////////////////////////////////////////////////////////////
	//  The user_function_compiler
	///////////////////////////////////////////////////////////////////////////
	bool user_function_compiler::operator()(double x) const {
		program.save_code(op_data);
		program.save_data(x);
		return true;
	}

	bool user_function_compiler::operator()(ast::variable const& x) const {
		program.save_code(op_load);
		program.save_variable_name(x.var);
		return true;
	}

	bool user_function_compiler::operator()(ast::vector const& x) const {
		size_t n = x.values.size();
		std::shared_ptr<Matrix<double>> new_matrix;
		if (n != 0) {
			std::list<double>::const_iterator val = x.values.begin();
			new_matrix = std::shared_ptr<Matrix<double>>(new Matrix<double>(1, n));
			for (int column = 0; val != x.values.end(); ++val, ++column) {
				(*new_matrix)[0][column] = *val;
			}
		}
		else {
			new_matrix = std::shared_ptr<Matrix<double>>(new Matrix<double>());
		}
		program.save_code(op_data);
		program.save_data(new_matrix);
		return true;
	}

	bool user_function_compiler::operator()(ast::matrix const& x) const {
		size_t m = x.rhs.size();
		std::shared_ptr<Matrix<double>> new_matrix;
		if (m != 0) {
			std::list<client::ast::vector>::const_iterator vec = x.rhs.begin();
			size_t n = vec->values.size();
			new_matrix = std::shared_ptr<Matrix<double>>(new Matrix<double>(m, n));
			for (int row = 0; vec != x.rhs.end(); ++vec, ++row) {
				if (vec->values.size() != n) {
					std::cout << "矩阵列数不一致" << std::endl;
					return false;
				}
				std::list<double>::const_iterator val = vec->values.begin();
				for (int column = 0; val != vec->values.end(); ++val, ++column) {
					(*new_matrix)[row][column] = *val;
				}
			}
		}
		else {
			new_matrix = std::shared_ptr<Matrix<double>>(new Matrix<double>());
		}
		program.save_code(op_data);
		program.save_data(new_matrix);
		return true;
	}

	bool user_function_compiler::operator()(ast::operation const& x) const {
		if (!boost::apply_visitor(*this, x.operand_))
			return false;
		switch (x.operator_) {
		case '+': program.save_code(op_add); break;
		case '-': program.save_code(op_sub); break;
		case '*': program.save_code(op_mul); break;
		case '/': program.save_code(op_div); break;
		case '^': program.save_code(op_pow); break;
		default: BOOST_ASSERT(0); return false;
		}
		return true;
	}

	bool user_function_compiler::operator()(ast::signed_ const& x) const {
		if (!boost::apply_visitor(*this, x.operand_))
			return false;
		switch (x.sign) {
		case '-': program.save_code(op_neg); break;
		case '+': break;
		default: BOOST_ASSERT(0); return false;
		}
		return true;
	}

	bool user_function_compiler::operator()(ast::function_call const& x)const {
		program.save_code(op_func_call);
		program.save_variable_name(x.function_name);
		BOOST_FOREACH(ast::parameter_type const& param, x.args) {
			program.save_code(op_param);
			if (!boost::apply_visitor(*this, param)) {
				return false;
			}
		}
		return true;
	}

	bool user_function_compiler::operator()(ast::expression const& x) const {
		if (!boost::apply_visitor(*this, x.first))
			return false;
		BOOST_FOREACH(ast::operation const& oper, x.rest) {
			if (!(*this)(oper))
				return false;
		}
		return true;
	}

	bool user_function_compiler::operator()(ast::assignment const& x) const {
		if (!(*this)(x.rhs))return false;
		program.save_code(op_store);
		program.save_variable_name(x.lhs);
		return true;
	}

	bool user_function_compiler::operator()(ast::variable_declaration const& x) const {
		if (x.rhs.size()>0) {//存在赋值语句
			BOOST_FOREACH(ast::expression const& expr, x.rhs) {
				if (!(*this)(expr))
					return false;
			}
			for (auto &elem : x.vars) {
				program.save_code(op_store_var);
				program.save_variable_name(elem);
			}
		}
		else {
			for (auto &elem : x.vars) {
				program.declara_var(elem);
			}
		}
		return true;
	}

	bool user_function_compiler::operator()(ast::define_function const& x) const {

		return true;
	}

	bool user_function_compiler::operator()(ast::statement_list const& x) const {
		BOOST_FOREACH(ast::statement const& statement, x) {
			if (!boost::apply_visitor(*this, statement)) {
				return false;
			}
		}
		return true;
	}
}

