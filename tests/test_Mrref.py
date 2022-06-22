# -*- coding: utf-8 -*-
import unittest
import numpy as np
from algebric_functions_computation_Mrref import compute_M_rref, rref

class VerificationMrrefCalculationsAndFixedFreeVariablesSolutions(unittest.TestCase):
    def testCalculationMrref(self):
        p = 1
        order = 3
        bool_constraints = False
        save_Mrref = False
        save_pos_free_variables = False
        
        A = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        B = np.array([-1, 3, 5])
        M = compute_M_rref(A, B, bool_constraints,order, p, save_Mrref, save_pos_free_variables)
        expectedMatrix = [[1, 0, -1, 0], [0, 1, 2, 0], [0, 0, 0, 1]]
        
        for i in range(len(expectedMatrix)):
            for j in range(len(expectedMatrix[0])):
                assert(np.abs(M[i][j] - expectedMatrix[i][j]) < 10**-12)
        
        A = np.array([[1,1,1,1,1], [1,1,0,0,0], [0,1,1,0,0], [1,2,1,0,0], [2,2,1,1,1]])
        B = np.array([1,2,3,5,3])
        M = compute_M_rref(A, B, bool_constraints,order, p, save_Mrref, save_pos_free_variables)
        expectedMatrix = [[1, 0, 0, 1, 1, -2], [0, 1, 0, -1, -1, 4], [0, 0, 1, 1, 1, -1], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]
        for i in range(len(expectedMatrix)):
            for j in range(len(expectedMatrix[0])):
                assert(np.abs(M[i][j] - expectedMatrix[i][j]) < 10**-12)
            
        """ compute_M_rref is working but not rref, it is weird    
        B = np.array(B)
        B = B.reshape(len(B), 1)
        M = np.hstack([A,B])
        M, _, _ = rref(M)
        print(M)
        """
if __name__ == '__main__':
    unittest.main()         