import sympy
from honeycomb import HoneyCombMesh
from expressions_functions import cleanExpr
from geometrical_elements import hexCycle
from WSF_to_generate_regular_hexagon import WSFTogenerateRegularHexagon

def getExpr(p, order,
            f2proj
            ):
    """
    The purpose here is to project a func (f2proj) in the Wachspress basis 
    for the regular hexagon
    """
    # pitch = distance between two parallel sides

    # create hexagonal mesh with 1 hexagonal cell
    mesh = HoneyCombMesh(1, p)

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
            funcDict['w%s' % (i+1)] = func(i+1, p)[0](x,y)
    elif order == 2:
        coeffs = []
        for i in range(6):
            wi, wiip1,coeffs_1 = func(i+1, p)
            funcDict['w%s' % (i+1)] = wi(x,y)
            nextI = wsf1.hexCycle(i+2)
            funcDict['w%s%s' % (i+1, nextI)] = wiip1(x,y)
            coeffs+=coeffs_1
        unknowns = coeffs
        print(unknowns)
    elif order == 3:
        coeffs= []
        coeffs_i = []
        coeffs_ip1 = []
        unknowns = []
        
        for i in range(6):
            #avec third_order2
            wi,  wiiip1, wip1ip1i, coeffs_2,coeffs_iiip1,coeffs_ip1ip1i = func(i+1, p)
            
            #avec third_order_1
            #wi,  wiiip1, wip1ip1i = func(i+1, p, coordsDict, linesDict)
            
            funcDict['w%s' % (i+1)] = wi(x,y)
            nextI = hexCycle(i+2)
            funcDict['w%s%s%s' % (nextI, nextI, i+1)] = wip1ip1i(x,y)
            funcDict['w%s%s%s' % (i+1, i+1, nextI)] = wiiip1(x,y)
            
            #'''pour third order2
            #We'll store all the symbolic coefficients in a precise order
            coeffs = coeffs + coeffs_2
            coeffs_i = coeffs_i + coeffs_iiip1
            coeffs_ip1 = coeffs_ip1 + coeffs_ip1ip1i
            #'''  
        #pour third_order2
        unknowns = coeffs + coeffs_i + coeffs_ip1
       
    elif order == 4:
        #We'll get all the symbolic coefficients in one list called unknowns
        coeffs= []
        coeffs_i = []
        coeffs_ip1 = []
        coeffs_ip1ip1 = []
        unknowns = []
        
        for i in range(6):
            wi,wiiip1,wiip1,wip1ip1i,coeffs_3,coeffs_iiip1,coeffs_iip1, coeffs_ip1ip1i = func(i+1, p)
            funcDict['w%s'% (i+1)] = wi(x,y)
            nextI = hexCycle(i+2)
            funcDict['w%s%s%s'%(i+1,i+1,nextI)] = wiiip1(x,y)
            funcDict['w%s%s'% (i+1,nextI)] = wiip1(x,y)
            funcDict['w%s%s%s'%(nextI,nextI,i+1)] = wip1ip1i(x,y)
            
            #We'll store all the symbolic coefficients in a precise order
            coeffs = coeffs + coeffs_3
            coeffs_i = coeffs_i + coeffs_iiip1
            coeffs_ip1 = coeffs_ip1 + coeffs_iip1
            coeffs_ip1ip1 = coeffs_ip1ip1 + coeffs_ip1ip1i
            
        unknowns = coeffs + coeffs_i + coeffs_ip1 + coeffs_ip1ip1
    elif order == 5:
        #We'll get all the symbolic coefficients in one list called unknowns
        coeffs_4 = []
        coeffs_3_1 = []
        coeffs_3_2 = []
        coeffs_3_3 = []
        coeffs_3_4 = []
        unknowns = []
        
        for i in range(6):
            wi,wi_1,wi_2,wi_3,wi_4, coeffs_i,coeffs_i_1,coeffs_i_2,coeffs_i_3, coeffs_i_4  = func(i+1, p)
            funcDict['w%s'% (i+1)] = wi(x,y)
            nextI = hexCycle(i+2)
            funcDict['w%s_1'%(i+1)] = wi_1(x,y)
            funcDict['w%s_2'%(i+1)] = wi_2(x,y)
            funcDict['w%s_3'%(i+1)] = wi_3(x,y)
            funcDict['w%s_4'%(i+1)] = wi_4(x,y)
            
            #We'll store all the symbolic coefficients in a precise order
            coeffs_4 = coeffs_4 + coeffs_i
            coeffs_3_1 = coeffs_3_1 + coeffs_i_1
            coeffs_3_2 = coeffs_3_2 + coeffs_i_2
            coeffs_3_3 = coeffs_3_3 + coeffs_i_3
            coeffs_3_4 = coeffs_3_4 + coeffs_i_4
        unknowns = coeffs_4 + coeffs_3_1 + coeffs_3_2 + coeffs_3_3 + coeffs_3_4
        
    else:
        raise ValueError("Order not implemented")

    # For order 1, there is only one basis function per i and wi[0] returns that function
    # For orders 2 and 3, there will be 2 and 3 functions
    # corresponding to the various basis functions as described by
    # Gout and the Physor2020 paper, cf WachspressShapeFunctions class
    # with computeSecondOrderFunction and computeThirdOrderFunction
    # for more details
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
            nextI = wsf1.hexCycle(i+2)
            sumyiWi += cleanExpr(funcDict['w%s%s' % (i+1, nextI)] * f2proj(*coordsDict['a%s%s' % (i+1, nextI)]))
    elif order == 3:
        print(funcDict['w112'])
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

    if(order == 5 or order == 4 or order == 3 or order == 2):
        return sumyiWi, unknowns
    else:
        return sumyiWi



