import sympy
from irregular_pentagon import IrregularPentagon
from expressions_functions import cleanExpr
from WSF_to_generate_irregular_pentagon import WSFTogenerateIrregularPentagon
from geometrical_elements import pentaCycle

def getExpr(order,
            f2proj):
    """
    The purpose here is to project a func (f2proj) in the Wachspress basis 
    for the irregular pentagon
    """
    # pitch = distance between two parallel sides
    # p = 1.

    # create hexagonal mesh with 1 hexagonal cell
    mesh = IrregularPentagon()

    # order = 3
    # Build the WSF object to generate shape functions
    wsf1 = WSFTogenerateIrregularPentagon(order, mesh)

    # Choose the appropriate method to build the shape functions for a
    # given order
    if order == 3:
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
    coordsDict, linesDict = wsf1.coords, wsf1.lines
    
    if order == 3:
        coeffs= []
        coeffs_i = []
        coeffs_ip1 = []
        unknowns = []
        
        for i in range(5):
            #avec third_order2
            wi,  wiiip1, wip1ip1i, coeffs_2,coeffs_iiip1,coeffs_ip1ip1i = func(i)
            
            #avec third_order_1
            #wi,  wiiip1, wip1ip1i = func(i+1, p, coordsDict, linesDict)
            
            funcDict['w%s' % (i)] = wi(x,y)
            nextI = pentaCycle(i+1)
            funcDict['w%s%s%s' % (nextI, nextI, i)] = wip1ip1i(x,y)
            funcDict['w%s%s%s' % (i, i, nextI)] = wiiip1(x,y)
            
            #'''pour third order2
            #We'll store all the symbolic coefficients in a precise order
            coeffs = coeffs + coeffs_2
            coeffs_i = coeffs_i + coeffs_iiip1
            coeffs_ip1 = coeffs_ip1 + coeffs_ip1ip1i
            #'''  
        unknowns = coeffs + coeffs_i + coeffs_ip1     
    elif order == 4:
        #We'll get all the symbolic coefficients in one list called unknowns
        coeffs= []
        coeffs_i = []
        coeffs_ip1 = []
        coeffs_ip1ip1 = []
        unknowns = []
        
        for i in range(5):
            wi,wiiip1,wiip1,wip1ip1i,coeffs_3,coeffs_iiip1,coeffs_iip1, coeffs_ip1ip1i = func(i)
            funcDict['w%s'% (i)] = wi(x,y)
            nextI = pentaCycle(i+1)
            funcDict['w%s%s%s'%(i,i,nextI)] = wiiip1(x,y)
            funcDict['w%s%s'% (i,nextI)] = wiip1(x,y)
            funcDict['w%s%s%s'%(nextI,nextI,i)] = wip1ip1i(x,y)
            
            #We'll store all the symbolic coefficients in a precise order
            coeffs = coeffs + coeffs_3
            coeffs_i = coeffs_i + coeffs_iiip1
            coeffs_ip1 = coeffs_ip1 + coeffs_iip1
            coeffs_ip1ip1 = coeffs_ip1ip1 + coeffs_ip1ip1i
            
        unknowns = coeffs + coeffs_i + coeffs_ip1 + coeffs_ip1ip1
    else:
        raise ValueError("Order not implemented")
    
    if order == 3:
        print(funcDict['w001'])
        sumyiWi = cleanExpr(funcDict['w0'] * f2proj(*coordsDict['a0']))
        sumyiWi += cleanExpr(funcDict['w110'] * f2proj(*coordsDict['a110']))
        sumyiWi += cleanExpr(funcDict['w001'] * f2proj(*coordsDict['a001']))
        for i in range(1,5):
            sumyiWi += cleanExpr(funcDict['w%s' % (i)] * f2proj(*coordsDict['a%s' % (i)]))
            nextI = pentaCycle(i+1)
            sumyiWi += cleanExpr(funcDict['w%s%s%s' % (nextI, nextI, i)] * f2proj(*coordsDict['a%s%s%s' % (nextI, nextI, i)]))
            sumyiWi += cleanExpr(funcDict['w%s%s%s' % (i, i, nextI)] * f2proj(*coordsDict['a%s%s%s' % (i, i, nextI)]))
    elif order == 4:
        sumyiWi = cleanExpr(funcDict['w0'] * f2proj(*coordsDict['a0'])) + cleanExpr(funcDict['w001'] * f2proj(*coordsDict['a001'])) + cleanExpr(funcDict['w01'] * f2proj(*coordsDict['a01'])) + cleanExpr(funcDict['w110'] * f2proj(*coordsDict['a110']))
        for i in range(1,5):
            nextI = pentaCycle(i+1)
            sumyiWi += cleanExpr(funcDict['w%s' % (i)] * f2proj(*coordsDict['a%s' % (i)])) + cleanExpr(funcDict['w%s%s%s'%(i,i,nextI)] * f2proj(*coordsDict['a%s%s%s'%(i,i,nextI)])) + cleanExpr(funcDict['w%s%s'%(i,nextI)] * f2proj(*coordsDict['a%s%s'%(i,nextI)])) + cleanExpr(funcDict['w%s%s%s'%(nextI,nextI,i)] * f2proj(*coordsDict['a%s%s%s'%(nextI,nextI,i)]))
    else:
        raise ValueError("Order {} not implemented".format(order))

    return sumyiWi,unknowns
