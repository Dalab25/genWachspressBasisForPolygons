# -*- coding: utf-8 -*-

import sympy
import numpy as np
from honeycomb import HoneyCombMesh
from WSF_to_generate_regular_hexagon import WSFTogenerateRegularHexagon
from expressions_functions_regular_hexagon import getExpr
from expressions_functions import exprDictToExprDictDict, putAllNumeratorsOnCommonDenominator, getDictFromExpr, getDenominatorsForXiYj
from geometrical_elements import hexCycle, computeLi, computeConic, computeCubic, computeQuartic
from algebric_functions import setUpLinearSystemFromDictNew
from save_files_functions import save_solution_txt, save_matrices

def get_coeffs_Gout(order,p):
    
    
    mesh = HoneyCombMesh(1, p)

    # order = 3
    # Build the WSF object to generate shape functions
    wsf1 = WSFTogenerateRegularHexagon(order, mesh)

    # Get coordinates and geometrical objects
    _, linesDict = wsf1._getCoordsAndLines()
    coeffs_gout = wsf1.get_coeffs_Gout(linesDict,p)
    return coeffs_gout  
        

def getWSF(p, order):
    """
    Function to get the shape functions from the WSF basis
    """
    # pitch = distance between two parallel sides
    # p = 1.

    # create hexagonal mesh with 1 hexagonal cell
    mesh = HoneyCombMesh(1, p)

    # order = 3
    # Build the WSF object to generate shape functions
    wsf1 = WSFTogenerateRegularHexagon(order, mesh)

    # Get coordinates and geometrical objects
    coordsDict, linesDict = wsf1.getCoordsAndLines()

    # Choose the appropriate method to build the shape functions for a
    # given order
    if order == 1:
        func = wsf1.computeFirstOrderFunction
    elif order == 2:
        func = wsf1.computeSecondOrderFunction
    elif order == 3:
        func = wsf1.computeThirdOrderFunction
    elif order == 4:
        func = wsf1.computeFourthOrderFunction
    elif order == 5:
        func = wsf1.computeFifthOrderFunction
    elif order == 6:
        func = wsf1.computeSixthOrderFunction
    else:
        raise ValueError("Order %s is not yet available" % order)

    return func, coordsDict, linesDict

def order3WithAllSetsCoeffs(pitch, bool_constraints = False,bool_save_matrices = False,bool_save_solution = False):
    """
    This function contains the main procedure to compute the missing coefficients 
    for the r^i, r^i_j functions at order 3.
    The idea is to project every monomial that belongs to the WSF basis
    To then create a linear system where the unknowns are the missing coefficients
    Hence, the linear system is solved and so the coefficients are obtained by the least-square method at first.
    """
    x = sympy.Symbol('x')
    y = sympy.Symbol('y')

    order = 3
    dictAllMonomialsProjected = {}
    
    degrees = [(i, j) for i in range(order+1) for j in range(order+1) if i+j <= order]
    for beta, gamma in degrees:
        print("-----> Checking beta = {}, gamma = {} ==> {}*{}".format(beta, gamma, x**beta, y**gamma))
        func = lambda a, b: a**beta * b**gamma
        
        expr,unknowns = getExpr(pitch, order, func)
        
        func, coords,_ = getWSF(pitch, order)
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
    
        if(beta == 0 and gamma == 0):
            for j in range(1,7):
                nextJ = hexCycle(j+1)
                
                ajjjp1x, ajjjp1y = coords['a%d%d%d' % (j,j,nextJ)][0], coords['a%d%d%d' % (j,j,nextJ)][1]
                ajp1jp1jx, ajp1jp1jy = coords['a%d%d%d' % (nextJ,nextJ,j)][0], coords['a%d%d%d' % (nextJ,nextJ,j)][1]
            
                if(j==1):
                    aj_1j_1jx, aj_1j_1jy = coords['a%d%d%d' % (6,6,j)][0], coords['a%d%d%d' % (6,6,j)][1]
                    ajjj_1x, ajjj_1y = coords['a%d%d%d' % (j,j,6)][0], coords['a%d%d%d' % (j,j,6)][1]
                else:    
                    aj_1j_1jx, aj_1j_1jy = coords['a%d%d%d' % (j-1,j-1,j)][0], coords['a%d%d%d' % (j-1,j-1,j)][1]
                    ajjj_1x, ajjj_1y = coords['a%d%d%d' % (j,j,j-1)][0], coords['a%d%d%d' % (j,j,j-1)][1]
                
                #coefficients pour les contraintes pour la fonction conique
                c3jjjp1 = computeConic(unknowns[6*(j-1):6*j], ajjjp1x, ajjjp1y)
                c3jp1jp1j = computeConic(unknowns[6*(j-1):6*j], ajp1jp1jx, ajp1jp1jy)
                c3j_1j_1j = computeConic(unknowns[6*(j-1):6*j], aj_1j_1jx, aj_1j_1jy)
                c3jjj_1 = computeConic(unknowns[6*(j-1):6*j], ajjj_1x, ajjj_1y)
                #print("unknows balayage 10: ", unknowns[10*(j-1):10*j])
                
                #coefficients pour les contraintes la droite pour les points d'ancrage a_iiip1
                c2j_1_jp1jp1j = computeLi(unknowns[36 + (j-1)*3:36 + (j*3)], ajp1jp1jx, ajp1jp1jy)
                #print("unknows balayage 60: ", unknowns[60 + 6*(j-1):60 + 6*j])
                
                #coefficients pour les contraintes de la droite les points d'ancrage a_ip1ip1i
                c2j_3_jjjp1 = computeLi(unknowns[54 + (j-1)*3:54 + (j*3)], ajjjp1x, ajjjp1y)
                #print("unknows balayage 96: ", unknowns[96 + 6*(j-1):96 + 6*j])
            
                #print("c_i ", c3jjjp1)
                #print("c_2", c3jjp1)
                #petit test sur lk
                
                lk = "supp%d%d%d"%(j,nextJ,j)
                keys = list(dictCoeff.keys())
                if lk in keys:
                    raise ValueError("Blup")
                
                if bool_constraints == True:
                    dictCoeff["supp%d%d%d"%(j,nextJ,j)] = c3jjjp1
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+6)] = c3jp1jp1j
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+12)] = c3j_1j_1j
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+18)] = c3jjj_1
                    
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+24)] = c2j_1_jp1jp1j
                    
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+30)] = c2j_3_jjjp1

        
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

    A, B = setUpLinearSystemFromDictNew(dictAllMonomialsProjected, unknowns, pitch)
    solution, res, _, _ = np.linalg.lstsq(A, B)
    print(">>> unknowns = \n", unknowns)
    print("solution: ",solution)
    print("residual: ",res)
    #'''
    if(bool_save_matrices == True):
        save_matrices(A,B,order,pitch,6, bool_constraints)
    if(bool_save_solution == True):
        save_solution_txt(bool_constraints, order, pitch,6, solution)
    
    return solution

def order4WithAllSetsCoeffs(pitch, bool_constraints = False,bool_save_matrices = False, bool_save_solution = False):
    """
    This function contains the main procedure to compute the missing coefficients 
    for the r^i, r^i_j functions at order 4.
    The idea is to project every monomial that belongs to the WSF basis
    To then create a linear system where the unknowns are the missing coefficients
    Hence, the linear system is solved and so the coefficients are obtained by the least-square method at first.
    """
    x = sympy.Symbol('x')
    y = sympy.Symbol('y')

    order = 4

    dictAllMonomialsProjected = {}
    
    degrees = [(i, j) for i in range(order+1) for j in range(order+1) if i+j <= order]
    for beta, gamma in degrees:
        print("-----> Checking beta = {}, gamma = {} ==> {}*{}".format(beta, gamma, x**beta, y**gamma))
        func = lambda a, b: a**beta * b**gamma
        
        expr,unknowns = getExpr(pitch, order, func)
        
        func, coords, _ = getWSF(pitch, order)
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

        if(beta == 0 and gamma == 0):
            for j in range(1,7):
                nextJ = hexCycle(j+1)
                
                ajjjp1x, ajjjp1y = coords['a%d%d%d' % (j,j,nextJ)][0], coords['a%d%d%d' % (j,j,nextJ)][1]
                ajjp1x, ajjp1y = coords['a%d%d' % (j,nextJ)][0], coords['a%d%d' % (j,nextJ)][1]
                ajp1jp1jx, ajp1jp1jy = coords['a%d%d%d' % (nextJ,nextJ,j)][0], coords['a%d%d%d' % (nextJ,nextJ,j)][1]
            
                if(j==1):
                    aj_1j_1jx, aj_1j_1jy = coords['a%d%d%d' % (6,6,j)][0], coords['a%d%d%d' % (6,6,j)][1]
                    aj_1jx, aj_1jy = coords['a%d%d' % (6,j)][0], coords['a%d%d' % (6,j)][1]
                    ajjj_1x, ajjj_1y = coords['a%d%d%d' % (j,j,6)][0], coords['a%d%d%d' % (j,j,6)][1]
                else:    
                    aj_1j_1jx, aj_1j_1jy = coords['a%d%d%d' % (j-1,j-1,j)][0], coords['a%d%d%d' % (j-1,j-1,j)][1]
                    aj_1jx, aj_1jy = coords['a%d%d' % (j-1,j)][0], coords['a%d%d' % (j-1,j)][1]
                    ajjj_1x, ajjj_1y = coords['a%d%d%d' % (j,j,j-1)][0], coords['a%d%d%d' % (j,j,j-1)][1]

            
                #coefficients pour les contraintes pour la fonction cubique
                c3jjjp1 = computeCubic(unknowns[10*(j-1):10*j], ajjjp1x, ajjjp1y)
                c3jjp1 = computeCubic(unknowns[10*(j-1):10*j], ajjp1x, ajjp1y)
                c3jp1jp1j = computeCubic(unknowns[10*(j-1):10*j], ajp1jp1jx, ajp1jp1jy)
                c3j_1j_1j = computeCubic(unknowns[10*(j-1):10*j], aj_1j_1jx, aj_1j_1jy)
                c3j_1j = computeCubic(unknowns[10*(j-1):10*j], aj_1jx, aj_1jy)
                c3jjj_1 = computeCubic(unknowns[10*(j-1):10*j], ajjj_1x, ajjj_1y)
                #print("unknows balayage 10: ", unknowns[10*(j-1):10*j])
                
                #coefficients pour les contraintes de la fonction conique pour les points d'ancrage a_iiip1
                c2j_1_jjp1 = computeConic(unknowns[60 + (j-1)*6:60 + (j*6)], ajjp1x, ajjp1y)
                c2j_1_jp1jp1j = computeConic(unknowns[60 + (j-1)*6:60 + (j*6)], ajp1jp1jx, ajp1jp1jy)
                #print("unknows balayage 60: ", unknowns[60 + 6*(j-1):60 + 6*j])
                
                #coefficients pour les contraintes de la fonction conique pour les points d'ancrage a_iiip1
                c2j_2_jjjp1 = computeConic(unknowns[96 + (j-1)*6:96 + (j*6)], ajjjp1x, ajjjp1y)
                c2j_2_jp1jp1j = computeConic(unknowns[96 + (j-1)*6:96 + (j*6)], ajp1jp1jx, ajp1jp1jy)
                #print("unknows balayage 60: ", unknowns[60 + 6*(j-1):60 + 6*j])
                
                
                #coefficients pour les contraintes de la fonction conique pour les points d'ancrage a_ip1ip1i
                c2j_3_jjjp1 = computeConic(unknowns[132 + (j-1)*6:132 + (j*6)], ajjjp1x, ajjjp1y)
                c2j_3_jjp1 = computeConic(unknowns[132 + (j-1)*6:132 + (j*6)], ajjp1x, ajjp1y)
                #print("unknows balayage 96: ", unknowns[96 + 6*(j-1):96 + 6*j])
            
                #print("c_i ", c3jjjp1)
                #print("c_2", c3jjp1)
                #petit test sur lk
                
                lk = "supp%d%d%d"%(j,nextJ,j)
                keys = list(dictCoeff.keys())
                if lk in keys:
                    raise ValueError("Blup")
                
                if bool_constraints == True:
                    dictCoeff["supp%d%d%d"%(j,nextJ,j)] = c3jjjp1
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+6)] = c3jjp1
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+12)] = c3jp1jp1j
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+18)] = c3j_1j_1j
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+24)] = c3j_1j
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+30)] = c3jjj_1
                    
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+36)] = c2j_1_jjp1
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+42)] = c2j_1_jp1jp1j
                    
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+48)] = c2j_2_jjjp1
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+54)] = c2j_2_jp1jp1j
                    
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+60)] = c2j_3_jjjp1
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+66)] = c2j_3_jjp1
            
        
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
    A, B = setUpLinearSystemFromDictNew(dictAllMonomialsProjected, unknowns, pitch)
    solution,res,rank,S = np.linalg.lstsq(A, B)
    print(">>> unknowns = \n", unknowns)
    print("solution: ",solution)
    print("residual: ",res)
    #'''
    
    if(bool_save_matrices == True):
        save_matrices(A,B,order,pitch,6, bool_constraints)
    if(bool_save_solution == True):
        save_solution_txt(bool_constraints, order, pitch,6, solution)
       
    return solution

def order5WithAllSetsCoeffs(pitch, bool_constraints = False,bool_save_matrices = True, bool_save_solution = True):
    """
    This function contains the main procedure to compute the missing coefficients 
    for the r^i, r^i_j functions at order 5.
    The idea is to project every monomial that belongs to the WSF basis
    To then create a linear system where the unknowns are the missing coefficients
    Hence, the linear system is solved and so the coefficients are obtained by the least-square method at first.
    """
    x = sympy.Symbol('x')
    y = sympy.Symbol('y')

    order = 5

    dictAllMonomialsProjected = {}
    
    #test de resolution pour y**4
    #beta = 0
    #gamma = 4
    degrees = [(i, j) for i in range(order+1) for j in range(order+1) if i+j <= order]
    for beta, gamma in degrees:
        print("-----> Checking beta = {}, gamma = {} ==> {}*{}".format(beta, gamma, x**beta, y**gamma))
        func = lambda a, b: a**beta * b**gamma
        
        expr,unknowns = getExpr(pitch, order, func)
        
        func, coords, lines = getWSF(pitch, order)
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
    
        if(bool_constraints == True):
            for j in range(1,7):
                nextJ = hexCycle(j+1)
                
                #coeffs_i = unknowns[(j-1)*order*(order+1)/2:j*order*(order+1)/2]
                coeffs_i = unknowns[(j-1)*15:j*15]
                ajx, ajy= coords['a%d' % j][0], coords['a%d' %j][1]
                aj1_x, aj1_y = coords['a%d1' %j][0], coords['a%d1' %j][1]
                aj2_x, aj2_y = coords['a%d2' %j][0], coords['a%d2' %j][1]
                aj3_x, aj3_y = coords['a%d3' %j][0], coords['a%d3' %j][1]
                aj4_x, aj4_y = coords['a%d4' %j][0], coords['a%d4' %j][1]

                if(j==1):
                    aj_1_1_x, aj_1_1_y = coords['a61'][0], coords['a61'][1]
                    aj_1_2_x, aj_1_2_y = coords['a62'][0], coords['a62'][1]
                    aj_1_3_x, aj_1_3_y = coords['a63'][0], coords['a63'][1]
                    aj_1_4_x, aj_1_4_y = coords['a64'][0], coords['a64'][1]
                    
                else:    
                    j_1 = j - 1
                    aj_1_1_x, aj_1_1_y = coords['a%d1' %j_1][0], coords['a%d1' %j_1][1]
                    aj_1_2_x, aj_1_2_y = coords['a%d2' %j_1][0], coords['a%d2' %j_1][1]
                    aj_1_3_x, aj_1_3_y = coords['a%d3' %j_1][0], coords['a%d3' %j_1][1]
                    aj_1_4_x, aj_1_4_y = coords['a%d4' %j_1][0], coords['a%d4' %j_1][1]
                    
                #coefficients pour les contraintes pour la fonction quadratique
                c4j_1 = computeQuartic(coeffs_i, aj1_x, aj1_y)
                c4j_2 = computeQuartic(coeffs_i, aj2_x, aj2_y)
                c4j_3 = computeQuartic(coeffs_i, aj3_x, aj3_y)
                c4j_4 = computeQuartic(coeffs_i, aj4_x, aj4_y)
                c4j_1_1 = computeQuartic(coeffs_i, aj_1_1_x, aj_1_1_y)
                c4j_1_2 = computeQuartic(coeffs_i, aj_1_2_x, aj_1_2_y)
                c4j_1_3 = computeQuartic(coeffs_i, aj_1_3_x, aj_1_3_y)
                c4j_1_4 = computeQuartic(coeffs_i, aj_1_4_x, aj_1_4_y)
                
                #s = 6*order*(order+1)/2
                #coeffs_i_1 = unknowns[s + (j-1)*order*(order-1)/2:s + j*order*(order-1)/2]
                coeffs_i_1 = unknowns[90 + 10*(j-1): 90 + 10*j]
                #coefficients pour les contraintes de la fonction cubique pour les points d'ancrage a_j1
                c3j1_2 = computeCubic(coeffs_i_1, aj2_x, aj2_y)
                c3j1_3 = computeCubic(coeffs_i_1, aj3_x, aj3_y)
                c3j1_4 = computeCubic(coeffs_i_1, aj4_x, aj4_y)

                #s = s + 6*order*(order-1)/2
                #coeffs_i_2 = unknowns[s + (j-1)*order*(order-1)/2:s + j*order*(order-1)/2]
                coeffs_i_2 = unknowns[150 + 10*(j-1): 150 + 10*j]
                #coefficients pour les contraintes de la fonction cubique pour les points d'ancrage a_j2
                c3j2_1 = computeCubic(coeffs_i_2, aj1_x, aj1_y)
                c3j2_3 = computeCubic(coeffs_i_2, aj3_x, aj3_y)
                c3j2_4 = computeCubic(coeffs_i_2, aj4_x, aj4_y)
                
                #s = s + 6*order*(order-1)/2
                #coeffs_i_3 = unknowns[s + (j-1)*order*(order-1)/2:s + j*order*(order-1)/2]
                coeffs_i_3 = unknowns[210 + 10*(j-1): 210 + 10*j]
                #coefficients pour les contraintes de la fonction cubique pour les points d'ancrage a_j3
                c3j3_1 = computeCubic(coeffs_i_3, aj1_x, aj1_y)
                c3j3_2 = computeCubic(coeffs_i_3, aj2_x, aj2_y)
                c3j3_4 = computeCubic(coeffs_i_3, aj4_x, aj4_y)
                
                #s = s + 6*order*(order-1)/2
                #coeffs_i_4 = unknowns[s + (j-1)*order*(order-1)/2:s + j*order*(order-1)/2]
                coeffs_i_4 = unknowns[270 + 10*(j-1): 270 + 10*j]
                #coefficients pour les contraintes de la fonction cubique pour les points d'ancrage a_j3
                c3j4_1 = computeCubic(coeffs_i_4, aj1_x, aj1_y)
                c3j4_2 = computeCubic(coeffs_i_4, aj2_x, aj2_y)
                c3j4_3 = computeCubic(coeffs_i_4, aj3_x, aj3_y)

                
                lk = "supp%d%d%d"%(j,nextJ,j)
                keys = list(dictCoeff.keys())
                if lk in keys:
                    raise ValueError("Blup")
                
                if bool_constraints == True:
                    dictCoeff["supp%d%d%d"%(j,nextJ,j)] = c4j_1
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+6)] = c4j_2
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+12)] = c4j_3
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+18)] = c4j_4
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+24)] = c4j_1_1
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+30)] = c4j_1_2
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+36)] = c4j_1_3
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+42)] = c4j_1_4
                    
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+48)] = c3j1_2
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+54)] = c3j1_3
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+60)] = c3j1_4
                    
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+66)] = c3j2_1
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+72)] = c3j2_3
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+78)] = c3j2_4
                    
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+84)] = c3j3_1
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+90)] = c3j3_2
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+96)] = c3j3_4
                    
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+102)] = c3j4_1
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+108)] = c3j4_2
                    dictCoeff["supp%d%d%d"%(j,nextJ,j+114)] = c3j4_3
        
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
    A, B = setUpLinearSystemFromDictNew(dictAllMonomialsProjected, unknowns, pitch)
    solution,res,rank,S = np.linalg.lstsq(A, B)
    print(">>> unknowns = \n", unknowns)
    print("solution: ",solution)
    print("residual: ",res)
    #'''
    
    if(bool_save_matrices == True):
        save_matrices(A,B,order,pitch,6, bool_constraints)
    if(bool_save_solution == True):
        save_solution_txt(bool_constraints, order, pitch,6, solution)
        
    return solution

def main():

    bool_constraints = False
    bool_save_matrices = True
    bool_save_solution = True
    p_ = [1]
    
    for p in p_:
        order3WithAllSetsCoeffs(p, bool_constraints,bool_save_matrices,bool_save_solution)
        order4WithAllSetsCoeffs(p, bool_constraints,bool_save_matrices, bool_save_solution)
        order5WithAllSetsCoeffs(p, bool_constraints,bool_save_matrices,bool_save_solution)

    print("The coefficients for the WSF are generated")


if __name__ == "__main__":
    main()
