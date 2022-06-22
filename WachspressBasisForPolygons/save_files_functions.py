import numpy as np
import os

def save_matrices(A, B, order, pitch, nSides, constraints):
    """
    Function to save the left and right members from the linear system to get the missing coefficients
    It is useful to create the rref matrix
    """
    os.makedirs('../results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/matrices/matrice_A', exist_ok=True)
    os.makedirs('../results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/matrices/matrice_B', exist_ok=True)
 
    if constraints == True:
        save_matrice_A = '../results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/matrices/matrice_A/matrice_A_constraints' + str(pitch)
        save_matrice_B = '../results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/matrices/matrice_B/matrice_B_constraints' + str(pitch)
    else:
        print("Saving")
        save_matrice_A = '../results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/matrices/matrice_A/matrice_A_no_constraints' + str(pitch)
        save_matrice_B = '../results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/matrices/matrice_B/matrice_B_no_constraints' + str(pitch)
   
    np.savetxt(save_matrice_A + '.txt',A)
    np.savetxt(save_matrice_B + '.txt',B)
    
def save_solution_txt(constraints, order, pitch, nSides, solution):
    """
    Function to save the computed coefficients for the r^i, r^i_j functions
    """
    os.makedirs('../results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/symbolic_coefficients/', exist_ok=True)
    if constraints == False:
        np.savetxt('../results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/symbolic_coefficients/solution_order_' + str(order) + '_no_constraints' + str(pitch) + '.txt', solution)
    else:   
        np.savetxt('../results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/symbolic_coefficients/solution_order_' + str(order) + '_constraints' + str(pitch) + '.txt', solution)
    
def save_M_rref_txt(order,pitch, M_rref, constraints=False, nSides=6):
    """
    Function to save the rref matrix
    """
    os.makedirs('../results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/matrices/M_rref', exist_ok=True)
    if constraints == False:
        save_M_rref = '../results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/matrices/M_rref/M_rref_order_' + str(order) + '_p_' + str(pitch) + '_no_constraints'
    else:
        save_M_rref = '../results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/matrices/M_rref/M_rref_order_' + str(order) + '_p_' + str(pitch) + '_constraints'
    np.savetxt(save_M_rref + '.txt', M_rref)
  
def save_pos_free_variables(order, pitch, pos_free_variables, nSides=6):
    """
    Function to save the placement of the free variables
    """
    os.makedirs('../results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/pos_free_variables', exist_ok=True)
    np.savetxt('../results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/pos_free_variables/pos_free_variables_order_' + str(order) + '_p_' + str(pitch) + '.txt', pos_free_variables)
    
    
    
    