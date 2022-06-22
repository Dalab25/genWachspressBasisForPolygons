# -*- coding: utf-8 -*-
import numpy as np
from algebric_functions_computation_Mrref import compute_M_rref, computationSolWithFreeVariables0
from convert_data_into_variables import convert_data_into_matrix_A_B


def main():
    p = 1
    order = 3
    nSides = 6
    method_rref = 0
    bool_constraints = False
    save_Mrref = True
    save_pos_free_variables = True
    #Computation of the rref matrix
    A,B,_ = convert_data_into_matrix_A_B(order,p, bool_constraints, nSides)
    M_rref1 = compute_M_rref(A, B, bool_constraints,order, p, save_Mrref, save_pos_free_variables, method_rref, nSides)
    #print(M_rref1, np.linalg.matrix_rank(M_rref1))
    
    #Computing the solution with the free variables
    computationSolWithFreeVariables0(order, p, bool_constraints, nSides)
    
    #Test to ensure that the linear systems build with and without the geometrical constraints are equivalent
    #recreate_sol_x_d_with_x_l_c_or_x_l_sc(order, p, bool_constraints, compare_solution = False)
    
if __name__ == "__main__":
    main()
    
