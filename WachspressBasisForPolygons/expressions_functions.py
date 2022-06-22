import sympy
import numpy as np

def expr2Dict(expr, *x, verbose=0):
    """
    From an expression, the idea is to create a dict whose keys
    are monomials. The values are sympy expression that came from 
    the projection of the monomial inside the WSF basis
    """
    collectedMonomials = sympy.Poly(expr, *x).as_expr()
    if verbose == 1:
        print(collectedMonomials)
    returndict = {}
    for i in sympy.Add.make_args(collectedMonomials):
        if verbose > 1:
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            print(i)
        term, coeff = i.as_independent(*x)[::-1]
        if verbose > 1:
            print(term, coeff)
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        keys = list(returndict)
        if term in keys:
            raise ValueError("bam")
            returndict[term].append(coeff)
        else:
            coeffs = []
            coeffs.append(coeff)
            returndict[term] = coeffs
            # for k, v in returndict.items():
            #     print(k, v)
    # print(returndict)
    return returndict


def expr2DictNoList(expr, *x, verbose=0):
    collectedMonomials = sympy.Poly(expr, *x).as_expr()
    if verbose == 1:
        print(collectedMonomials)
    returndict = {}
    for i in sympy.Add.make_args(collectedMonomials):
        if verbose > 1:
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            print(i)
        term, coeff = i.as_independent(*x)[::-1]
        if verbose > 1:
            print(term, coeff)
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        keys = list(returndict)
        if term in keys:
            raise ValueError(">>>> key in dict", term)
        else:
            returndict[term] = coeff
            if verbose > 1:
                for k, v in returndict.items():
                    print(k, v)
        # input(">")
    # print(returndict)
    return returndict


def cleanExpr(expr):
    """
    Function to factorize the sympy expression
    """
    # val = expr
    val = sympy.expand(expr)
    # Define symbolic x and y
    # x = sympy.Symbol('x')
    # y = sympy.Symbol('y')
    # val = sympy.Poly(expr, x, y)
    return val


def getDictFromExpr(expr, x, y, verbose=0):
    """
    Returns dictCoeff from a given expression.
    """
    dictCoeff = {}
    colly = sympy.collect(expr, y)
    # print(colly)
    for i in range(10):
        yi = colly.coeff(y, i)
        collx = sympy.collect(yi, x)
        for j in range(10):
            xj = collx.coeff(x, j)
            if xj == 0:
                pass
            else:
                dictCoeff[y**i*x**j] = xj
                # print(i, j, xj)

    if verbose:
        for k, v in dictCoeff.items():
            print(k, v)

    return dictCoeff

def getDenominatorsForXiYj(dictCoeff, verbose=0):
    """
    Returns the denominators in same order as keys read from dictCoeff or provided.
    """
    alldens = {}
    allnums = {}
    for key, value in dictCoeff.items():
        if verbose:
            print(key, "\t", value, value.args)
        if isinstance(value, sympy.Float):
            monomials = [value]
        else:
            #input(value.args)
            monomials = value.args
        nums = []
        dens = []
        alldensi = []
        for val in monomials:
            num, den = sympy.fraction(val)
            if verbose > 1:
                print("num = ", num, "\t den = ",den)
            nums.append(num)
            if den not in dens:
                dens.append(den)
                alldensi.append(den)
        alldens[key] = alldensi
        allnums[key] = nums

    if verbose:
        for i, (k,d) in enumerate(alldens.items()):
            print(i, k, d)

    return alldens


def putAllNumeratorsOnCommonDenominator(dictCoeff, alldens, verbose):
    dictCoeff2 = {}
    for ik, (key, value) in enumerate(dictCoeff.items()):
        if isinstance(value, sympy.Float):
            monomials = [value]
        else:
            monomials = value.args
        numsForMono = []
        for val in monomials:
            num, den = sympy.fraction(val)
            if verbose > 1:
                print(">>> Start common den for:")
                print(num, den)
                for dik in alldens[key]:
                    print(">", dik)
            lcm = np.lcm.reduce(alldens[key])
            numsForMono.append(num*lcm)
        dictCoeff2[key] = sum(numsForMono)

    if verbose:
        for k, v in dictCoeff2.items():
            print(k, v)

    return dictCoeff2


def exprDictToExprDictDict(dictCoeff,
                           unknowns,
                           fixup=True,
                           verbose=0
                           ):
    """
    Transforms a dictionary with keys x**n*y**m and values being
    expressions of the unknown coefficients to a dict of dict of the
    form:
    dict[x**n*y**m] = { a_i: ..., b_i: ...., c_i: ... }
    """
    newDict = {}
    for key, value in dictCoeff.items():
        if verbose:
            print(">>>>> Processing: ", key, value)
        num, den = sympy.fraction(value)
        #print("ATTENTION", *unknowns)
        unkCoeffDict = expr2DictNoList(num, *unknowns)
        unknowsInEq = unkCoeffDict.keys()
        notInKeys = set(unknowns) - set(unknowsInEq)
        zeros = [0] * len(notInKeys)
        notInKeysDict = dict(zip(notInKeys, zeros))
        unkCoeffDict.update(notInKeysDict)
        newDict[key] = unkCoeffDict

    if fixup:
        epsilon = 1.e-12
        for i, (key, equation) in enumerate(newDict.items()):
            # print(">>>>> key", key, equation)
            for j, u in enumerate(unknowns):
                coeff = equation[u]
                if abs(coeff) < epsilon:
                    equation[u] = 0
                coeff = equation[u]
            try:
                coeff = equation[1]
                if abs(coeff) < epsilon:
                    equation[1] = 0
                coeff = equation[1]
            except:
                pass

    if verbose:
        for i, (key, equation) in enumerate(newDict.items()):
            print(">>>>> key", key, equation)

    return newDict

