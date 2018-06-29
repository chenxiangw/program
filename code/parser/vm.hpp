
#pragma once

#include"compiler.h"
#include <vector>

namespace client{
	enum byte_code{
        op_neg, 
        op_add,   
        op_sub,  
        op_mul,  
        op_div,   
		op_pow,

        op_load_var,    //  load a variable
		op_load_mat,    //  load a matrix
		op_store,
        op_store_var,   //  store a variable
		op_store_var_plural, //  store many variables
		op_store_mat,//  store a matrix
        op_dou,     
		op_mat,

		op_param_val,
		op_param_mat,
		op_func_call
    };

    ///////////////////////////////////////////////////////////////////////////
    //  The Virtual Machine
    ///////////////////////////////////////////////////////////////////////////
    
    class vmachine{

    public:
        vmachine(std::map<std::string, std::pair<bool, double>> &variables,
		std::map < std::string, std::shared_ptr<Matrix<double>>> &matrixes,
		std::map < std::string, void(*)(client::vmachine  *vm)> &builtin_functions)
        : matrixes(matrixes)
        , variables(variables)
		, builtin_functions(builtin_functions)
        {
        }

		bool execute(client::code_gen::program &program);
		std::vector<
			boost::variant<double, std::shared_ptr<Matrix<double>>>
		> & get_data() { return data; }

		void clear() { data.clear(); }
		void print_stack();

    private:
		std::vector<boost::variant<double,
			std::shared_ptr<Matrix<double>>
		>> data;
		std::map<std::string, std::pair<bool, double>> &variables;
		std::map < std::string, std::shared_ptr<Matrix<double>>> &matrixes;
		std::map < std::string, void(*)(client::vmachine  *vm) > &builtin_functions;
    };

}


