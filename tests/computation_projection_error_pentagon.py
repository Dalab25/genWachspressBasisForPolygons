import sympy
import numpy as np
from geometrical_elements import computeAdjointForIrregularPentagon, pentaCycle
from expressions_functions import getDictFromExpr, cleanExpr
from irregular_pentagon import IrregularPentagon
from pentagonal_functions_generalised_basis import PentagonalFunctionsGeneralisedBasis


def getExprToComputeTheProjectionError(pentagon, order, f2proj, boolConstraints):
    """
    This function multiplies every Wachspress functions in order to identify
    x**beta y**gamma *adjoint(x,y) = WSF(x,y)
    """
    p = 1
    # pitch = distance between two parallel sides
    # p = 1.
    if(order > 2):
        if(boolConstraints):
            coeffsFile = "../results/polygon_5sides/order" + str(order) + "/symbolic_coefficients/solution_order_" + str(order) + "_constraints1.txt"
        else:
            coeffsFile = "../results/polygon_5sides/order" + str(order) + "/symbolic_coefficients/solution_order_" + str(order) + "_no_constraints1.txt"   

    else:
        coeffsFile = ""
    wsf1 = PentagonalFunctionsGeneralisedBasis(order,  pentagon,  p, coeffsFile)

    # Get coordinates and geometrical objects
    coordsDict, _ = wsf1.coords, wsf1.lines
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
    else:
        raise ValueError("Order %s is not yet available" % order)

    # Define symbolic x and y
    x = sympy.Symbol('x')
    y = sympy.Symbol('y')

    # Get the shape functions.
    # func(i+1, p, coordsDict, linesDict) returns a tuple with 1, 2 or
    # 3 basis functions depending on the order
    funcDict = {}
    
    if order == 1:
        for i in range(5):
            funcDict['w%s' % (i)] = func(i)(x,y)#*computeAdjointForIrregularPentagon(wsf1.pentagon.coeffsAdjoint, x, y)
   
    elif order == 2:
        for i in range(5):
            wi, wiip1 = func(i)
            funcDict['w%s' % (i)] = wi(x,y)#*computeAdjointForIrregularPentagon(wsf1.pentagon.coeffsAdjoint, x, y)
   
            nextI = pentaCycle(i+1)
            funcDict['w%s%s' % (i, nextI)] = wiip1(x,y)#*computeAdjointForIrregularPentagon(wsf1.pentagon.coeffsAdjoint, x, y)
   
    
    elif order == 3:
        for i in range(5):
            wi, wiiip1, wip1ip1i = func(i)
            funcDict['w%s' % (i)] = wi(x,y)#*computeAdjointForIrregularPentagon(wsf1.pentagon.coeffsAdjoint, x, y))
            nextI = pentaCycle(i+1)
            funcDict['w%s%s%s' % (i, i, nextI)] = wiiip1(x,y)#*computeAdjointForIrregularPentagon(wsf1.pentagon.coeffsAdjoint, x, y))
            funcDict['w%s%s%s' % (nextI, nextI, i)] = wip1ip1i(x,y)#*computeAdjointForIrregularPentagon(wsf1.pentagon.coeffsAdjoint, x, y))
            
            
    elif order == 4:
        for i in range(5):
            wi,wiiip1,wiip1,wip1ip1i = func(i)
            funcDict['w%s'% (i)] = wi(x,y)#*computeAdjointForIrregularPentagon(wsf1.pentagon.coeffsAdjoint, x, y)
            nextI = pentaCycle(i+1)
            funcDict['w%s%s%s'%(i,i,nextI)] = wiiip1(x,y)#*computeAdjointForIrregularPentagon(wsf1.pentagon.coeffsAdjoint, x, y)
            funcDict['w%s%s'% (i,nextI)] = wiip1(x,y)#*computeAdjointForIrregularPentagon(wsf1.pentagon.coeffsAdjoint, x, y)
            funcDict['w%s%s%s'%(nextI,nextI,i)] = wip1ip1i(x,y)#*computeAdjointForIrregularPentagon(wsf1.pentagon.coeffsAdjoint, x, y)
            
    else:
        raise ValueError("Order not implemented")
    
    sumyiWi = 0
    if order == 1:
        for i in range(5):
            sumyiWi += cleanExpr(funcDict['w%s' % (i)] * f2proj(*coordsDict['a%s' % (i)]))
    elif order == 2:
        for i in range(5):
            sumyiWi += cleanExpr(funcDict['w%s' % (i)] * f2proj(*coordsDict['a%s' % (i)]))
            nextI = pentaCycle(i+1)
            sumyiWi += cleanExpr(funcDict['w%s%s' % (i, nextI)] * f2proj(*coordsDict['a%s%s' % (i, nextI)]))
    elif order == 3:
        for i in range(5):
            sumyiWi += cleanExpr(funcDict['w%s' % (i)] * f2proj(*coordsDict['a%s' % (i)]))
            nextI = pentaCycle(i+1)
            sumyiWi += cleanExpr(funcDict['w%s%s%s' % (nextI, nextI, i)] * f2proj(*coordsDict['a%s%s%s' % (nextI, nextI, i)]))
            sumyiWi += cleanExpr(funcDict['w%s%s%s' % (i, i, nextI)] * f2proj(*coordsDict['a%s%s%s' % (i, i, nextI)]))
    
    elif order == 4:
        for i in range(5):
            nextI = pentaCycle(i+1)
            sumyiWi += cleanExpr(funcDict['w%s' % (i)] * f2proj(*coordsDict['a%s' % (i)])) + cleanExpr(funcDict['w%s%s%s'%(i,i,nextI)] * f2proj(*coordsDict['a%s%s%s'%(i,i,nextI)])) + cleanExpr(funcDict['w%s%s'%(i,nextI)] * f2proj(*coordsDict['a%s%s'%(i,nextI)])) + cleanExpr(funcDict['w%s%s%s'%(nextI,nextI,i)] * f2proj(*coordsDict['a%s%s%s'%(nextI,nextI,i)]))
    else:
        raise ValueError("Order {} not implemented".format(order))

    return sumyiWi


def calcul_err(order, boolConstraints):
    ''' General method 
    #we calculate through the basis the projection on a monomial
    #example with the projection on x*y
    
    #manual approximation
    appro_x = 1.0000001*x
    appro_y = 1.0000001*y
    
    #calculation of the errors
    err_ = x*y - appro_x*appro_y
    
    #get the coefficient in front of the monomials trough a dictionnary
    dict = getDictFromExpr(err_,x,y)
    
    #list of the coefficients that we'll normalize
    L_coeff_err = []
    for i in dict:
        print(dict[i]) #print the coefficient in front of the dictionnary
        L_coeff_err.append(dict[i])
    
    #we get the error through the normalization of this list
    print(np.linalg.norm(L_coeff_err))
    '''
    x = sympy.Symbol("x")
    y = sympy.Symbol("y")
    L_err = []
    pentagon = IrregularPentagon()
    
    degrees = [(i, j) for i in range(order+1) for j in range(order+1) if i+j <= order]

    for gamma, beta in degrees:
        func = lambda a, b: (a**gamma) * (b**beta)
        sumwi = getExprToComputeTheProjectionError(pentagon, order,func,boolConstraints)
        
        adjoint = computeAdjointForIrregularPentagon(pentagon.coeffsAdjoint,x,y)
        func = func(x,y)*adjoint #multiplying by the adjoint
        func = cleanExpr(func)
            
        #
        expr_err = func - sumwi

        #Creating a dictionnary to identify for every monomials
        dictCoeff = getDictFromExpr(expr_err, x, y)
        
        #Creating a list for the projection of the monomials number i
        L_err_i = []
        for i in dictCoeff:
            print(dictCoeff[i])
            L_err_i.append(float(dictCoeff[i]))
        #print(L_err_i)
    
        #Adding the error of the projection number i in the list
        L_err.append(np.linalg.norm(L_err_i))
        
    print(L_err)
    print("Error : ", np.linalg.norm(L_err)) 
        

def main():
    boolConstraints = True
    order = 3
    p = 1
    calcul_err(order, boolConstraints)

    """
    x = sympy.Symbol("x")
    y = sympy.Symbol("y")
    
    # pitch = distance between two parallel sides
    # p = 1.
    if(order > 2):
        if(boolConstraints):
            coeffsFile = "../results/polygon_5sides/order" + str(order) + "/symbolic_coefficients/solution_order_" + str(order) + "_constraints1.txt"
        else:
            coeffsFile = "../results/polygon_5sides/order" + str(order) + "/symbolic_coefficients/solution_order_" + str(order) + "_no_constraints1.txt"   

    else:
        coeffsFile = ""
    """
    """
    pentagon = IrregularPentagon()
    wsf1 = PentagonalFunctionsGeneralisedBasis(order,  pentagon,  p, coeffsFile)
    (wi,wiiip1, wip1ip1i) = wsf1.computeThirdOrderFunction(0)
    exp = (wi(x,y)*computeAdjointForIrregularPentagon(wsf1.pentagon.coeffsAdjoint, x, y)).evalf()
    print("Expression", exp)
     """
if __name__ == "__main__":
    main()