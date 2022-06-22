import sympy
import numpy as np
from honeycomb import HoneyCombMesh
from hexagonal_functions_generalised_basis import HexagonalFunctionsGeneralisedBasis
from expressions_functions import cleanExpr
from triangular_partition_quadrature import TriangularPartitionQuadrature
from geometrical_elements import hexCycle

def getExprToComputeTheProjectionError(hexagon, order, p, f2proj, boolHomothetic):
    """
    This function multiplies every Wachspress functions in order to identify
    x**beta y**gamma *adjoint(x,y) = WSF(x,y)
    """
    # pitch = distance between two parallel sides
    # p = 1.
    if(order > 2):
        if(boolHomothetic):
            coeffsFile = "../results/polygon_6sides/order" + str(order) + "/symbolic_coefficients/solution_order_" + str(order) + "_pitch_1_free_variables_0.txt"
        else:
            coeffsFile = "../results/polygon_6sides/order" + str(order) + "/symbolic_coefficients/solution_order_" + str(order) + "_no_constraints" + str(p) + ".txt"
            #coeffsFile = "../results/polygon_6sides/order" + str(order) + "/symbolic_coefficients/solution_order_3_no_constraints32.90896534380867.txt"
        
    else:
        coeffsFile = ""
    wsf1 = HexagonalFunctionsGeneralisedBasis(hexagon, order, coeffsFile, homothetic=boolHomothetic)
    # Get coordinates and geometrical objects
    coordsDict, linesDict = wsf1.getCoordsLines()
    
    # Choose the appropriate method to build the shape functions for a given order
    if order == 1:
        func = wsf1.computeFirstOrderFunction
    elif order == 2:
        #func = wsf1.computeSecondOrderFunctionGout
        func = wsf1.computeSecondOrderFunction
    elif order == 3:
        func = wsf1.computeThirdOrderFunction
    elif order == 4:
        func = wsf1.computeFourthOrderFunction
    elif order == 5:
        func = wsf1.computeFifthOrderFunction
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
        for i in range(6):
            funcDict['w%s' % (i+1)] = func(i+1, p, coordsDict, linesDict)(x,y)
   
    elif order == 2:
        for i in range(6):
            wi, wiip1 = func(i+1, p, coordsDict, linesDict)
            funcDict['w%s' % (i+1)] = wi(x,y)
            nextI = hexCycle(i+2)
            funcDict['w%s%s' % (i+1, nextI)] = wiip1(x,y)
    
    elif order == 3:
        for i in range(6):
            wi,  wiiip1, wip1ip1i = func(i+1)
            funcDict['w%s' % (i+1)] = wi(x,y)
            nextI = hexCycle(i+2)
            funcDict['w%s%s%s' % (i+1, i+1, nextI)] = wiiip1(x,y)
            funcDict['w%s%s%s' % (nextI, nextI, i+1)] = wip1ip1i(x,y)
            
            '''
            for a_i in coordsDict:
                coordsDict[a_i] = np.array(coordsDict[a_i])
                #print("Avant", coordsDict[a_i][0])
                #print("Après", coordsDict[a_i][0]/(p**beta))
                coordsDict[a_i][0] /= (p**beta)
                coordsDict[a_i][1] /=(p**gamma) 
                coordsDict[a_i] = tuple(coordsDict[a_i])
            '''
            
    elif order == 4:
        for i in range(6):
            wi,wiiip1,wiip1,wip1ip1i = func(i+1)
            funcDict['w%s'% (i+1)] = wi(x,y)
            nextI = hexCycle(i+2)
            funcDict['w%s%s%s'%(i+1,i+1,nextI)] = wiiip1(x,y)
            funcDict['w%s%s'% (i+1,nextI)] = wiip1(x,y)
            funcDict['w%s%s%s'%(nextI,nextI,i+1)] = wip1ip1i(x,y)
            
            '''
            for a_i in coordsDict:
                coordsDict[a_i] = np.array(coordsDict[a_i])
                #print("Avant", coordsDict[a_i][0])
                #print("Après", coordsDict[a_i][0]/(p**beta))
                coordsDict[a_i][0] /= (p**beta)
                coordsDict[a_i][1] /=(p**gamma) 
                coordsDict[a_i] = tuple(coordsDict[a_i])
            '''
    elif order == 5:        
        for i in range(6):
            wi,wi_1,wi_2,wi_3,wi_4  = func(i+1)
            funcDict['w%s'% (i+1)] = wi(x,y)
            nextI = hexCycle(i+2)
            funcDict['w%s_1'%(i+1)] = wi_1(x,y)
            funcDict['w%s_2'%(i+1)] = wi_2(x,y)
            funcDict['w%s_3'%(i+1)] = wi_3(x,y)
            funcDict['w%s_4'%(i+1)] = wi_4(x,y)
            
    else:
        raise ValueError("Order not implemented")

    order = len(funcDict.keys())//6
    if order == 1:
        sumyiWi = cleanExpr(funcDict['w1'] * f2proj(*coordsDict['a1']))
        for i in range(1,6):
            sumyiWi += cleanExpr(funcDict['w%s' % (i+1)] * f2proj(*coordsDict['a%s' % (i+1)]))
    elif order == 2:
        sumyiWi = cleanExpr(funcDict['w1'] * f2proj(*coordsDict['a1']))
        sumyiWi += cleanExpr(funcDict['w12'] * f2proj(*coordsDict['a12']))
        for i in range(1,6):
            sumyiWi += cleanExpr(funcDict['w%s' % (i+1)] * f2proj(*coordsDict['a%s' % (i+1)]))
            nextI = hexCycle(i+2)
            sumyiWi += cleanExpr(funcDict['w%s%s' % (i+1, nextI)] * f2proj(*coordsDict['a%s%s' % (i+1, nextI)]))
    elif order == 3:
        
        sumyiWi = cleanExpr(funcDict['w1'] * f2proj(*coordsDict['a1']))
        sumyiWi += cleanExpr(funcDict['w221'] * f2proj(*coordsDict['a221']))
        sumyiWi += cleanExpr(funcDict['w112'] * f2proj(*coordsDict['a112']))
        for i in range(1,6):
            sumyiWi += cleanExpr(funcDict['w%s' % (i+1)] * f2proj(*coordsDict['a%s' % (i+1)]))
            nextI = hexCycle(i+2)
            sumyiWi += cleanExpr(funcDict['w%s%s%s' % (nextI, nextI, i+1)] * f2proj(*coordsDict['a%s%s%s' % (nextI, nextI, i+1)]))
            sumyiWi += cleanExpr(funcDict['w%s%s%s' % (i+1, i+1, nextI)] * f2proj(*coordsDict['a%s%s%s' % (i+1, i+1, nextI)]))
    
    elif order == 4:
        sumyiWi = cleanExpr(funcDict['w1'] * f2proj(*coordsDict['a1'])) + cleanExpr(funcDict['w112'] * f2proj(*coordsDict['a112'])) + cleanExpr(funcDict['w12'] * f2proj(*coordsDict['a12'])) + cleanExpr(funcDict['w221'] * f2proj(*coordsDict['a221']))
        for i in range(1,6):
            nextI = hexCycle(i+2)
            sumyiWi += cleanExpr(funcDict['w%s' % (i+1)] * f2proj(*coordsDict['a%s' % (i+1)])) + cleanExpr(funcDict['w%s%s%s'%(i+1,i+1,nextI)] * f2proj(*coordsDict['a%s%s%s'%(i+1,i+1,nextI)])) + cleanExpr(funcDict['w%s%s'%(i+1,nextI)] * f2proj(*coordsDict['a%s%s'%(i+1,nextI)])) + cleanExpr(funcDict['w%s%s%s'%(nextI,nextI,i+1)] * f2proj(*coordsDict['a%s%s%s'%(nextI,nextI,i+1)]))
    elif order == 5:
        sumyiWi = cleanExpr(funcDict['w1'] * f2proj(*coordsDict['a1'])) + cleanExpr(funcDict['w1_1'] * f2proj(*coordsDict['a11'])) + cleanExpr(funcDict['w1_2'] * f2proj(*coordsDict['a12'])) + cleanExpr(funcDict['w1_3'] * f2proj(*coordsDict['a13'])) + cleanExpr(funcDict['w1_4'] * f2proj(*coordsDict['a14']))
        for i in range(1,6):
            nextI = hexCycle(i+2)
            sumyiWi += cleanExpr(funcDict['w%s' % (i+1)] * f2proj(*coordsDict['a%s' % (i+1)])) + cleanExpr(funcDict['w%s_1'%(i+1)] * f2proj(*coordsDict['a%s1'%(i+1)])) + cleanExpr(funcDict['w%s_2'%(i+1)] * f2proj(*coordsDict['a%s2'%(i+1)])) + cleanExpr(funcDict['w%s_3'%(i+1)] * f2proj(*coordsDict['a%s3'%(i+1)])) + cleanExpr(funcDict['w%s_4'%(i+1)] * f2proj(*coordsDict['a%s4'%(i+1)]))
    else:
        raise ValueError("Order {} not implemented".format(order))

    return sumyiWi

''' General method '''
'''   
    #we calculate through the base the projection on a monomial
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

def computationL2Error(order, p, vectorL2Err,boolHomothetic):
    """
    For every order:
        vecL2Error = []
        For every pitch:
            For all monomial x^gamma y^beta belonging to Pk:
                err_proj_L2 = integration(x^gamma y^beta - projWSF(x^gamma, y^beta))
                vecL2Error.append(sqrt(err_proj_L2))
        errL2OrderK = np.linalg.norm(vecL2Error)
        
    """
    hexagon = HoneyCombMesh(1,p)
    nQuad = 15
    quads = TriangularPartitionQuadrature.getNodesAndWeights(nQuad,p)
    degrees = [(i, j) for i in range(order+1) for j in range(order+1) if i+j <= order]
    x,y = sympy.Symbol("x"), sympy.Symbol("y")
    vector = []
    for gamma, beta in degrees:
        #Integration of every monomial on the hexagonal cell of pitch p
        err = 0
        func2proj = lambda a, b: (a**gamma) * (b**beta)
        sumwIxGammayBeta = getExprToComputeTheProjectionError(hexagon, order, p, func2proj,boolHomothetic)
        sumwIxGammayBeta = sympy.lambdify([x,y], sumwIxGammayBeta)
        for i in range(len(quads)):
            xi, yi, wi = quads[i][0], quads[i][1], quads[i][2]
            err += wi*(func2proj(xi,yi) - sumwIxGammayBeta(xi,yi))**2
        #Taking the L2 error of integration of the current monomial
        vector.append(np.sqrt(err))
    vectorL2Err.append(np.linalg.norm(err))
             
def main():
    orders = [3]
    p_ = [1]#[0.01, 0.1 1, 10, 100,1000,10000]
    boolHomothetic = True
    #orders = [3]
    #p_ = [1,10]
    for order in orders:
        vectorL2err = []
        for p in p_:
            computationL2Error(order, p, vectorL2err,boolHomothetic)
            print("L2 Error for order " + str(order) + ", for p =", p_)
            print(vectorL2err)
        print("L2 Norm for all pitch: ", np.linalg.norm(vectorL2err))
                



if __name__ == "__main__":
    main()