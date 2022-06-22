# -*- coding: utf-8 -*-

import sympy
import numpy as np
from irregular_pentagon import IrregularPentagon
from WSF_to_generate_irregular_pentagon import WSFTogenerateIrregularPentagon
from expressions_functions import exprDictToExprDictDict, putAllNumeratorsOnCommonDenominator, getDictFromExpr, getDenominatorsForXiYj
from expressions_functions_irregular_pentagon import getExpr
from algebric_functions import setUpLinearSystemFromDictNewIrregularPentagon
from save_files_functions import save_solution_txt, save_matrices
from geometrical_elements import pentaCycle, computeLi, computeConic, computeCubic

def getWSF(order):
    """
    Function to get the shape functions from the WSF basis
    """
    # pitch = distance between two parallel sides
    # p = 1.

    # create hexagonal mesh with 1 hexagonal cell
    mesh = IrregularPentagon()

    # order = 3
    # Build the WSF object to generate shape functions
    wsf1 = WSFTogenerateIrregularPentagon(order, mesh)

    # Get coordinates and geometrical objects
    coordsDict, linesDict = wsf1.coords, wsf1.lines
    
    if order == 3:
        func = wsf1.computeThirdOrderFunction
    elif order == 4:
        func = wsf1.computeFourthOrderFunction
    else:
        raise ValueError("Order %s is not yet available" % order)

    return func,coordsDict, linesDict



def order3WithAllSetsCoeffs(bool_constraints = False,bool_save_matrices = False,bool_save_solution = False):
    """
    This function contains the main procedure to compute the missing coefficients
    for the r^i, r^i_j functions at order 3.
    The idea is to project every monomial that belongs to the WSF basis
    To then create a linear system where the unknowns are the missing coefficients
    Hence, the linear system is solved and so the coefficients are obtained by the least-square method at first.
    """
    x = sympy.Symbol('x')
    y = sympy.Symbol('y')
    
    n = 5
    pitch = 1
    order = 3
    

    dictAllMonomialsProjected = {}
    
    degrees = [(i, j) for i in range(order+1) for j in range(order+1) if i+j <= order]
    for beta, gamma in degrees:
        print("-----> Checking beta = {}, gamma = {} ==> {}*{}".format(beta, gamma, x**beta, y**gamma))
        func = lambda a, b: a**beta * b**gamma
        
        expr,unknowns = getExpr(order, func)
        
        _, coords, _ = getWSF(order)
        expr = sympy.expand(expr)
        print("expr: ", expr)
            # for e in expr.args:
            #     print(e)
    
        dictCoeff = getDictFromExpr(expr, x, y)
        print("DictCoeff: ", dictCoeff)
        alldens = getDenominatorsForXiYj(dictCoeff)
        print("alldens:", alldens)
        dictXiYjTermsCommonDen = putAllNumeratorsOnCommonDenominator(dictCoeff,
                                                                         alldens,
                                                                         verbose=0
                                                                         )
    
        dictCoeff = dictXiYjTermsCommonDen
        print("DictCoeff2: ", dictCoeff)
    
        # for k, v in dictCoeff.items():
        #     print(k, v)
        if(bool_constraints):
            if(beta == 0 and gamma == 0):
                for j in range(n):
                    #nextJ = WachspressShapeFunctions.hexCycle(j+1)
                    nextJ = pentaCycle(j+1)
                    ajjjp1x, ajjjp1y = coords['a%d%d%d' % (j,j,nextJ)][0], coords['a%d%d%d' % (j,j,nextJ)][1]
                    ajp1jp1jx, ajp1jp1jy = coords['a%d%d%d' % (nextJ,nextJ,j)][0], coords['a%d%d%d' % (nextJ,nextJ,j)][1]
                
                    if(j==0):
                        aj_1j_1jx, aj_1j_1jy = coords['a%d%d%d' % (4,4,j)][0], coords['a%d%d%d' % (4,4,j)][1]
                        ajjj_1x, ajjj_1y = coords['a%d%d%d' % (j,j,4)][0], coords['a%d%d%d' % (j,j,4)][1]
                    else:    
                        h = j-1
                        aj_1j_1jx, aj_1j_1jy = coords['a%d%d%d' % (h,h,j)][0], coords['a%d%d%d' % (h,h,j)][1]
                        ajjj_1x, ajjj_1y = coords['a%d%d%d' % (j,j,h)][0], coords['a%d%d%d' % (j,j,h)][1]
                    
                    #Coefficients for the constraint on the conic function
                    c3jjjp1 = computeConic(unknowns[6*j:6*(j+1)], ajjjp1x, ajjjp1y)
                    c3jp1jp1j = computeConic(unknowns[6*j:6*(j+1)], ajp1jp1jx, ajp1jp1jy)
                    c3j_1j_1j = computeConic(unknowns[6*j:6*(j+1)], aj_1j_1jx, aj_1j_1jy)
                    c3jjj_1 = computeConic(unknowns[6*j:6*(j+1)], ajjj_1x, ajjj_1y)
    
                    #Coefficients for the line function anchored to the pint a_iiip1
                    c2j_1_jp1jp1j = computeLi(unknowns[n*6 + j*3:n*6 + 3*(j+1)], ajp1jp1jx, ajp1jp1jy)
    
                    #Coefficients for the line function anchored to the pint a_ip1ip1i
                    c2j_3_jjjp1 = computeLi(unknowns[n*6 + n*3 + 3*j:n*6 + n*3 + 3*(j+1)], ajjjp1x, ajjjp1y)
                    
                    
                    lk = "supp%d%d%d"%(j,nextJ,j)
                    keys = list(dictCoeff.keys())
                    if lk in keys:
                        raise ValueError("Blup")
                    
                        dictCoeff["supp%d%d%d"%(j,nextJ,j)] = c3jjjp1
                        dictCoeff["supp%d%d%d"%(j,nextJ,j+n)] = c3jp1jp1j
                        dictCoeff["supp%d%d%d"%(j,nextJ,j+2*n)] = c3j_1j_1j
                        dictCoeff["supp%d%d%d"%(j,nextJ,j+3*n)] = c3jjj_1
                        
                        dictCoeff["supp%d%d%d"%(j,nextJ,j+4*n)] = c2j_1_jp1jp1j
                        
                        dictCoeff["supp%d%d%d"%(j,nextJ,j+5*n)] = c2j_3_jjjp1
    
            
        newDict = exprDictToExprDictDict(dictCoeff, unknowns,
                                             fixup=True,
                                             verbose=0
                                             )
        dictAllMonomialsProjected[(beta, gamma)] = newDict

    for k, v in dictAllMonomialsProjected.items():
        print(k)
        for l, vv in v.items():
            print(l,vv)
            
    print("dictAllMonomialsProjected: ",dictAllMonomialsProjected)

    A, B = setUpLinearSystemFromDictNewIrregularPentagon(dictAllMonomialsProjected, unknowns)
    solution,res,rank,S = np.linalg.lstsq(A, B)
    print(">>> unknowns = \n", unknowns)
    print("solution: ",solution)
    print("residual: ",res)
    #'''
    if(bool_save_matrices == True):
        save_matrices(A, B, order, pitch, n, bool_constraints)
    if(bool_save_solution == True):
        save_solution_txt(bool_constraints, order, pitch, n, solution)
    
    return solution

def order4WithAllSetsCoeffs(bool_constraints = False,bool_save_matrices = False,bool_save_solution = False):
    """
    This function contains the main procedure to compute the missing coefficients 
    for the r^i, r^i_j functions at order 4.
    The idea is to project every monomial that belongs to the WSF basis
    To then create a linear system where the unknowns are the missing coefficients
    Hence, the linear system is solved and so the coefficients are obtained by the least-square method at first.
    """
    x = sympy.Symbol('x')
    y = sympy.Symbol('y')
    
    n = 5
    pitch = 1
    order = 4

    dictAllMonomialsProjected = {}
    
    degrees = [(i, j) for i in range(order+1) for j in range(order+1) if i+j <= order]
    for beta, gamma in degrees:
        print("-----> Checking beta = {}, gamma = {} ==> {}*{}".format(beta, gamma, x**beta, y**gamma))
        func = lambda a, b: a**beta * b**gamma
        
        expr,unknowns = getExpr(order, func)
        
        func, coords, _ = getWSF(order)
        expr = sympy.expand(expr)
        print("expr: ", expr)

        dictCoeff = getDictFromExpr(expr, x, y)
        print("DictCoeff: ", dictCoeff)
        alldens = getDenominatorsForXiYj(dictCoeff)
        print("alldens:", alldens)
        dictXiYjTermsCommonDen = putAllNumeratorsOnCommonDenominator(dictCoeff,
                                                                         alldens,
                                                                         verbose=0
                                                                         )
    
        dictCoeff = dictXiYjTermsCommonDen
        print("DictCoeff2: ", dictCoeff)
        
        if(bool_constraints):
            if(beta == 0 and gamma == 0):
                for j in range(5):
                    nextJ = pentaCycle(j+1)
                    
                    ajjjp1x, ajjjp1y = coords['a%d%d%d' % (j,j,nextJ)][0], coords['a%d%d%d' % (j,j,nextJ)][1]
                    ajjp1x, ajjp1y = coords['a%d%d' % (j,nextJ)][0], coords['a%d%d' % (j,nextJ)][1]
                    ajp1jp1jx, ajp1jp1jy = coords['a%d%d%d' % (nextJ,nextJ,j)][0], coords['a%d%d%d' % (nextJ,nextJ,j)][1]
                
                    if(j==0):
                        aj_1j_1jx, aj_1j_1jy = coords['a%d%d%d' % (4,4,j)][0], coords['a%d%d%d' % (4,4,j)][1]
                        aj_1jx, aj_1jy = coords['a%d%d' % (4,j)][0], coords['a%d%d' % (4,j)][1]
                        ajjj_1x, ajjj_1y = coords['a%d%d%d' % (j,j,4)][0], coords['a%d%d%d' % (j,j,4)][1]
                    else:    
                        h = j - 1
                        aj_1j_1jx, aj_1j_1jy = coords['a%d%d%d' % (h,h,j)][0], coords['a%d%d%d' % (h,h,j)][1]
                        aj_1jx, aj_1jy = coords['a%d%d' % (h,j)][0], coords['a%d%d' % (h,j)][1]
                        ajjj_1x, ajjj_1y = coords['a%d%d%d' % (j,j,h)][0], coords['a%d%d%d' % (j,j,h)][1]
    
    
                
                    #Coefficient for the cubic function
                    c3jjjp1 = computeCubic(unknowns[10*j:10*(j+1)], ajjjp1x, ajjjp1y)
                    c3jjp1 = computeCubic(unknowns[10*j:10*(j+1)], ajjp1x, ajjp1y)
                    c3jp1jp1j = computeCubic(unknowns[10*j:10*(j+1)], ajp1jp1jx, ajp1jp1jy)
                    c3j_1j_1j = computeCubic(unknowns[10*j:10*(j+1)], aj_1j_1jx, aj_1j_1jy)
                    c3j_1j = computeCubic(unknowns[10*j:10*(j+1)], aj_1jx, aj_1jy)
                    c3jjj_1 = computeCubic(unknowns[10*j:10*(j+1)], ajjj_1x, ajjj_1y)
                    #print("unknows balayage 10: ", unknowns[10*(j-1):10*j])
                    
                    #Coefficient for the conic function anchored to the pint a_iiip1
                    c2j_1_jjp1 = computeConic(unknowns[n*10 + j*6:n*10  + (j+1)*6], ajjp1x, ajjp1y)
                    c2j_1_jp1jp1j = computeConic(unknowns[n*10  + j*6: n*10 + (j+1)*6], ajp1jp1jx, ajp1jp1jy)
                    #print("unknows balayage 60: ", unknowns[60 + 6*(j-1):60 + 6*j])
                    
                    #Coefficients for the conic function anchored to the pint a_iip1
                    c2j_2_jjjp1 = computeConic(unknowns[n*10 + n*6 + j*6: n*10 + n*6 + (j+1)*6], ajjjp1x, ajjjp1y)
                    c2j_2_jp1jp1j = computeConic(unknowns[n*10 + n*6 + j*6: n*10 + n*6 + (j+1)*6], ajp1jp1jx, ajp1jp1jy)
                    #print("unknows balayage 60: ", unknowns[60 + 6*(j-1):60 + 6*j])
                    
                    #Coefficients for the conic function anchored to the pint a_ip1ip1i
                    c2j_3_jjjp1 = computeConic(unknowns[n*10 + 2*n*6  + j*6 : n*10 + 2*n*6+ (j+1)*6], ajjjp1x, ajjjp1y)
                    c2j_3_jjp1 = computeConic(unknowns[n*10 + 2*n*6 +  j*6 : n*10 + 2*n*6 + (j+1)*6], ajjp1x, ajjp1y)
    
                    lk = "supp%d%d%d"%(j,nextJ,j)
                    keys = list(dictCoeff.keys())
                    if lk in keys:
                        raise ValueError("Blup")
                    
                    if bool_constraints == True:
                        dictCoeff["supp%d%d%d"%(j,nextJ,j)] = c3jjjp1
                        dictCoeff["supp%d%d%d"%(j,nextJ,j+n)] = c3jjp1
                        dictCoeff["supp%d%d%d"%(j,nextJ,j+2*n)] = c3jp1jp1j
                        dictCoeff["supp%d%d%d"%(j,nextJ,j+3*n)] = c3j_1j_1j
                        dictCoeff["supp%d%d%d"%(j,nextJ,j+4*n)] = c3j_1j
                        dictCoeff["supp%d%d%d"%(j,nextJ,j+5*n)] = c3jjj_1
                        
                        dictCoeff["supp%d%d%d"%(j,nextJ,j+6*n)] = c2j_1_jjp1
                        dictCoeff["supp%d%d%d"%(j,nextJ,j+7*n)] = c2j_1_jp1jp1j
                        
                        dictCoeff["supp%d%d%d"%(j,nextJ,j+8*n)] = c2j_2_jjjp1
                        dictCoeff["supp%d%d%d"%(j,nextJ,j+9*n)] = c2j_2_jp1jp1j
                        
                        dictCoeff["supp%d%d%d"%(j,nextJ,j+10*n)] = c2j_3_jjjp1
                        dictCoeff["supp%d%d%d"%(j,nextJ,j+11*n)] = c2j_3_jjp1
            
        
        newDict = exprDictToExprDictDict(dictCoeff, unknowns,
                                             fixup=True,
                                             verbose=0
                                             )
        dictAllMonomialsProjected[(beta, gamma)] = newDict

    for k, v in dictAllMonomialsProjected.items():
        print(k)
        for l, vv in v.items():
            print(l,vv)
            
    print("dictAllMonomialsProjected: ",dictAllMonomialsProjected)
    #exit(">>> ok")
    #'''
    A, B = setUpLinearSystemFromDictNewIrregularPentagon(dictAllMonomialsProjected, unknowns)
    solution,res,rank,S = np.linalg.lstsq(A, B)
    print(">>> unknowns = \n", unknowns)
    print("solution: ",solution)
    print("residual: ",res)
    #'''
    
    if(bool_save_matrices == True):
        save_matrices(A, B, order, pitch, n, bool_constraints)
    if(bool_save_solution == True):
        save_solution_txt(bool_constraints, order, pitch, n, solution)
    
    return solution



def main():

    bool_constraints = False
    bool_save_matrices = True
    bool_save_solution = True
    order3WithAllSetsCoeffs(bool_constraints,bool_save_matrices,bool_save_solution)
    order4WithAllSetsCoeffs(bool_constraints,bool_save_matrices,bool_save_solution)


if __name__ == "__main__":
    main()
