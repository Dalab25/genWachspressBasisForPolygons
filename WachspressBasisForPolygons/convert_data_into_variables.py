import numpy as np

#nSides=6 by default because the regular hexagon is the main interest

def convert_data_into_matrix_A_B(order, pitch, bool_constraints,nSides=6):
    """
    Load the files A.txt and B.txt to turn them into matrix/vector
    """
    unknowns_i = (order*(order + 1)/2) + (order -1)*(order*(order-1))/2
    nb_cols = int(nSides*unknowns_i) 
    if(bool_constraints == True):
        A = np.loadtxt('../results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/matrices/matrice_A/matrice_A_constraints' + str(pitch) + '.txt', usecols = range(nb_cols))
        B = np.loadtxt('../results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/matrices/matrice_B/matrice_B_constraints' + str(pitch)+ '.txt')
        solution = np.loadtxt('../results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/symbolic_coefficients/solution_order_' + str(order) + '_constraints' + str(pitch) + '.txt')
    else:
        A = np.loadtxt('../results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/matrices/matrice_A/matrice_A_no_constraints' + str(pitch) + '.txt', usecols = range(nb_cols))
        B = np.loadtxt('../results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/matrices/matrice_B/matrice_B_no_constraints' + str(pitch)+ '.txt')
        solution = np.loadtxt('../results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/symbolic_coefficients/solution_order_' + str(order) + '_no_constraints' + str(pitch) + '.txt')
    return A, B,solution

def convert_data_into_matrix_M_rref(order, pitch, bool_constraints, nSides=6):
    """
    Load the file M_rref.txt to turn it into a matrix
    """
    unknowns_i = (order*(order + 1)/2) + (order -1)*(order*(order-1))/2
    nb_cols = int(nSides*unknowns_i) + 1
    if(bool_constraints == True):
        M_rref = np.loadtxt('../results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/matrices/M_rref/M_rref_order_' + str(order) + '_p_' + str(pitch) + '_constraints.txt', usecols = range(nb_cols))
    else:
        M_rref = np.loadtxt('../results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/matrices/M_rref/M_rref_order_' + str(order) + '_p_' + str(pitch) + '_no_constraints.txt', usecols = range(nb_cols))
    return M_rref

  
def convert_data_into_free_pos_variables(order, pitch, nSides=6):  
    """
    Load the file pos_free_variables.txt to turn it into a vector
    """
    pos_free_variables = np.loadtxt('../results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/pos_free_variables/pos_free_variables_order_' + str(order) + '_p_' + str(pitch) + '.txt',dtype ='float').astype(int)
    return pos_free_variables


