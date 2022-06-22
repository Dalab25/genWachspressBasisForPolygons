# -*- coding: utf-8 -*-

import sympy
import numpy as np
from itertools import combinations
from irregular_pentagon import IrregularPentagon

def backslash(A, b):
    """
    This is the equivalent of matlab A\b
    """
    num_vars = A.shape[1]
    rank = np.linalg.matrix_rank(A)
    if rank == num_vars:
        # print rank, num_vars
        # b = np.array(b)
        # r = np.linalg.matrix_rank(b)
        # print rank, r, A.ndim, b.ndim
        sol = np.linalg.lstsq(A, b, rcond=None)[0]  # not under-determined
        return sol
    else:
        it = 0
        for nz in combinations(range(num_vars), rank):
            it += 1
            # print(dir(combinations(range(num_vars), rank)))
            # print(range(num_vars))
            # exit()
            # the variables not set to zero
            try:
                sol = np.zeros((num_vars, 1))
                sol[nz, :] = np.asarray(np.linalg.solve(A[:, nz], b))
                return sol
            except np.linalg.LinAlgError:
                print("LinAlgError for rank = ", rank, " iter = ", it)
                pass  # picked bad variables, can't solve


def setUpLinearSystemFromDictNewIrregularPentagon(dictAllMono, unknowns):
    """
    Function to create a linear system from a given dictionnary
    The keys of the dictionnary are monomials
    The values are sympy expressions corresponding to the projection of the monomial in the WSF basis
    From there, an identification is done to create the linear system
    The identification is based on the expected values computed from the adjoint of the polygon
    """
    pentagon = IrregularPentagon()
    coeffsAdjoint = pentagon.coeffsAdjoint
    nEquationsList = []
    degrees = dictAllMono.keys()
    for k in degrees:
        nEquationsList.append(len(dictAllMono[k].keys()))
    

    nEquations = sum(nEquationsList)
    # nEquations = len(newDict.keys())  # nRows
    nUnknowns = len(unknowns)
    A = np.zeros(nEquations * nUnknowns).reshape(nEquations, nUnknowns)
    B = np.zeros(nEquations)
    A[:] = np.nan
    B[:] = np.nan

    x = sympy.Symbol('x')
    y = sympy.Symbol('y')

    for ideg, degree in enumerate(degrees):
        beta, gamma = degree
        print(">>> Identification for: ", x**(2+beta) * y**gamma, x**beta * y**(2+gamma), x**beta * y**gamma)

    for ideg, degree in enumerate(degrees):
        newDict = dictAllMono[degree]
        if ideg > 0:
            nPrevious = sum(nEquationsList[:ideg])
        else:
            nPrevious = 0
        beta, gamma = degree
        for i, (key, equation) in enumerate(newDict.items()):
            # print(beta, gamma, nPrevious, x**(2+beta) * y**gamma, x**beta * y**(2+gamma), x**beta * y**gamma, key, equation)
            if key in [x**(2+beta) * y**gamma]:
                B[nPrevious + i] = 1
            elif key in [x**beta * y**(2+gamma)]:
                B[nPrevious + i] = coeffsAdjoint[0]
            elif key in [x**(beta+1) * y**(gamma+1)]:
                B[nPrevious + i] = coeffsAdjoint[1]
            elif key in [x**(beta+1) * y**(gamma)]:
                B[nPrevious + i] = coeffsAdjoint[2]
            elif key in [x**(beta) * y**(gamma + 1)]:
                B[nPrevious + i] = coeffsAdjoint[3]
            elif key in [x**(beta) * y**(gamma)]:
                B[nPrevious + i] = coeffsAdjoint[4]
            else:
                B[nPrevious + i] = 0
            for j, u in enumerate(unknowns):
                coeff = equation[u]
                A[nPrevious + i,j] = coeff
            try:
                B[nPrevious + i] += -equation[1]
            except:
                pass

    print(">>> Matrix A = \n", A)
    print(">>> RHS B = \n", B)
    print(">>> Rank of A = ", np.linalg.matrix_rank(A))
    print(">>> nEquations = {}, nUnknowns = {}".format(nEquations, nUnknowns))

    return A, B


def setUpLinearSystemFromDictNew(dictAllMono, unknowns, pitch):
    """
    Function to create a linear system from a given dictionnary
    The keys of the dictionnary are monomials
    The values are sympy expressions corresponding to the projection of the monomial in the WSF basis
    From there, an identification is done to create the linear system
    The identification is based on the expected values computed from the adjoint of the polygon
    """
    nEquationsList = []
    degrees = dictAllMono.keys()
    for k in degrees:
        nEquationsList.append(len(dictAllMono[k].keys()))

    nEquations = sum(nEquationsList)
    # nEquations = len(newDict.keys())  # nRows
    nUnknowns = len(unknowns)
    A = np.zeros(nEquations * nUnknowns).reshape(nEquations, nUnknowns)
    B = np.zeros(nEquations)
    A[:] = np.nan
    B[:] = np.nan

    x = sympy.Symbol('x')
    y = sympy.Symbol('y')

    for ideg, degree in enumerate(degrees):
        beta, gamma = degree
        print(">>> Identification for: ", x**(2+beta) * y**gamma, x**beta * y**(2+gamma), x**beta * y**gamma)

    for ideg, degree in enumerate(degrees):
        newDict = dictAllMono[degree]
        if ideg > 0:
            nPrevious = sum(nEquationsList[:ideg])
        else:
            nPrevious = 0
        beta, gamma = degree
        for i, (key, equation) in enumerate(newDict.items()):
            # print(beta, gamma, nPrevious, x**(2+beta) * y**gamma, x**beta * y**(2+gamma), x**beta * y**gamma, key, equation)
            if key in [x**(2+beta) * y**gamma, x**beta * y**(2+gamma)]:
                B[nPrevious + i] = 1
            elif key in [x**beta * y**gamma]:
                B[nPrevious + i] = -pitch**2
            else:
                B[nPrevious + i] = 0
            for j, u in enumerate(unknowns):
                coeff = equation[u]
                A[nPrevious + i,j] = coeff
            try:
                B[nPrevious + i] += -equation[1]
            except:
                pass

    print(">>> Matrix A = \n", A)
    print(">>> RHS B = \n", B)
    print(">>> Rank of A = ", np.linalg.matrix_rank(A))
    print(">>> nEquations = {}, nUnknowns = {}".format(nEquations, nUnknowns))

    return A, B


def almostEqualRelativeAndAbs(a, b, maxDiff=1.e-12, maxRelDiff=1.e-11):
    """
    Check if numbers are close to each other according to a certain tolerance
    """
    # when comparing numbers near zero.
    diff = abs(a - b)
    if diff <= maxDiff:
        return True

    a = abs(a)
    b = abs(b)
    largest = b if b > a else a

    if diff <= largest * maxRelDiff:
        return True

    return False

def compareSolutionToRef(ref, sol, unknowns):
    for u in unknowns:
        refVal = ref[u]
        solVal = sol[u]
        print(u, refVal, solVal)
        assert(almostEqualRelativeAndAbs(refVal, solVal))

def calcul_residu(A,B,solution):
    res = B - np.dot(A,solution)
    res = np.linalg.norm(res)
    print(res)
