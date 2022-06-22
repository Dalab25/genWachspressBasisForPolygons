# -*- coding: utf-8 -*-
#
# File: wachspress_shape_functions.py
#
# from builtins import staticmethod
"""
This class defines the Wachspress shape functions and the geometric
elements required for their construction.
Reference: J. L. Gout, Rational Wachspress-type Finite Elements on
Regular Hexagons, Journal of Numerical Analysis (1985)
"""
import sympy
import numpy as np
import matplotlib.pyplot as plt
from hexagonal_nodal_dofs import buildCoordsAndLinesForWSF
from geometrical_elements import hexCycle, computeCircleForRegularHex, computeLi, computeConic, computeCubic, computeQuartic
from matplotlib.path import Path

"""
This class defines the Wachspress shape functions and the geometric
elements required for their construction.
To compute with the symbolic calulcation the basis, the denominator
is not taking into account.
"""

class WSFTogenerateRegularHexagon(object):

    def __init__(self, order, hexMesh):
        """ Constructor with a given order and hexMesh, required to
        obtain geometrical data such as lattice pitch, halfside of
        hexagon.  """
        self.__order = order
        self.__basisSize = 6 * self.__order
        self.__pitch = hexMesh.pitch
        self.__halfPitch = hexMesh.halfPitch
        self.__halfSide = hexMesh.halfSide
        self.__coords, self.__lines = buildCoordsAndLinesForWSF(self.__halfPitch, self.__halfSide, order)
        self.__func = None

        if self.__order == 1:
            self.__func = self.computeFirstOrderFunction
        elif self.__order == 2:
            self.__func = self.computeSecondOrderFunction
        elif self.__order == 3:
            self.__func = self.computeThirdOrderFunction
        elif self.__order == 4:
            self.__func = self.computeFourthOrderFunction
        elif self.__order == 5:
            self.__func = self.computeFifthOrderFunction
        else:
            raise ValueError("Order %s is not yet available" % order)

    @property
    def basisSize(self):
        """
        Returns the basis size for a given order.
        """
        return self.__basisSize

    def getCoordsAndLines(self):
        return self.__coords, self.__lines


    def computeFirstOrderFunction(self, j,p):
        """
        Function to build first-order Wachspress elements (cf. Gout 1985).
        """
        coords, lines = self.getCoordsAndLines()
        line = lambda i: lines['d%d' %   hexCycle(i)]
        func = lambda x, y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y)
        a0, a1 = coords['a%d' % j][0], coords['a%d' % j][1]
        qai = computeCircleForRegularHex(p, a0, a1)
        ci = qai / func(a0, a1)
        wi = lambda x, y: ci * func(x, y) / computeCircleForRegularHex(p, x, y)

        return (wi,)

    def computeSecondOrderFunction(self, j, p):
        """
        Function to build second-order Wachspress elements (cf. Gout 1985).
        """
        coords, lines = self.getCoordsAndLines()
        line = lambda i: lines['d%d' %   hexCycle(i)]
        linep = lines['dp%d' % j]
        func1 = lambda x, y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) * computeLi(linep, x, y)
        a0, a1 = coords['a%d' % j][0], coords['a%d' % j][1]
        qai = computeCircleForRegularHex(p, a0, a1)
        ci = qai / func1(a0, a1)
        wi = lambda x, y: ci * func1(x, y) / computeCircleForRegularHex(p, x, y)

        func2 = lambda x, y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) * computeLi(line(j), x, y)
        ip1 = lambda i:   hexCycle(i)
        aiip1 = coords['a%d%d' % (j, ip1(j+1))]
        qaiip1 = computeCircleForRegularHex(p, aiip1[0], aiip1[1])
        ciip1 = qaiip1 / func2(aiip1[0], aiip1[1])
        wiip1 = lambda x, y: ciip1 * func2(x, y) / computeCircleForRegularHex(p, x, y)

        (wi,wiip1)
    
    def computeThirdOrderFunction(self,j,p):
        coords,lines = self.getCoordsAndLines()
        line = lambda i: lines['d%d' %   hexCycle(i)]
        nextJ =   hexCycle(j+1)

        '''First terms of basis: wi'''

        #coefficients for the first function
        a_2 = sympy.Symbol('a_2_%d'%j)
        b_2 = sympy.Symbol('b_2_%d'%j)
        c_2 = sympy.Symbol('c_2_%d'%j)
        d_2 = sympy.Symbol('d_2_%d'%j)
        e_2 = sympy.Symbol('e_2_%d'%j)
        f_2 = sympy.Symbol('f_2_%d'%j)

        coeffs_2 = [a_2,b_2,c_2,d_2,e_2,f_2]

        func1 = lambda x,y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y)* computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) #* computeCubic(coeffs_3,x,y)
        func11 = lambda x, y: func1(x,y) * computeConic(coeffs_2, x, y)
        a0, a1 = coords['a%d' % j][0], coords['a%d' % j][1]
        qai = computeCircleForRegularHex(p, a0, a1)
        ci = qai / func1(a0, a1)
        wi = lambda x, y: func11(x, y) * ci #/ computeCircleForRegularHex(p, x, y) # * ci
        # wi = lambda x, y:  func1(x, y) * ci

        '''Second terms of basis: wiiip1'''

        #coefficients for the second function
        a_iiip1 = sympy.Symbol('a_iiip1_%d'%j)
        b_iiip1 = sympy.Symbol('b_iiip1_%d'%j)
        c_iiip1 = sympy.Symbol('c_iiip1_%d'%j)

        coeffs_iiip1 = [a_iiip1, b_iiip1, c_iiip1]
        func2 = lambda x,y:  computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) * computeLi(line(j), x, y) #* cls.computeConic(coeffs_iiip1, x,y)
        func22 = lambda x, y: func2(x,y) * computeLi(coeffs_iiip1, x, y)
        
        aiiip1 = coords['a%d%d%d' % (j,j,nextJ)]
        a0, a1 = aiiip1[0], aiiip1[1]
        qaiiip1 = computeCircleForRegularHex(p, a0, a1)
        ciiip1 = qaiiip1 / func2(a0, a1)
        wiiip1 = lambda x, y:  func22(x, y) * ciiip1 #/ computeCircleForRegularHex(p, x, y) # * ciiip1
        # wiiip1 = lambda x, y:  func2(x, y) * ciiip1

        '''Third terms of basis: wip1ip1i'''

        #Coefficients for the fourth function
        a_ip1ip1i = sympy.Symbol('a_ip1ip1i_%d'%j)
        b_ip1ip1i = sympy.Symbol('b_ip1ip1i_%d'%j)
        c_ip1ip1i = sympy.Symbol('c_ip1ip1i_%d'%j)


        coeffs_ip1ip1i = [a_ip1ip1i, b_ip1ip1i, c_ip1ip1i]
        func3 = lambda x,y:  computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) * computeLi(line(j), x, y) #* cls.computeConic(coeffs_ip1ip1i, x,y)
        func33 = lambda x, y: func3(x,y) * computeLi(coeffs_ip1ip1i, x, y)
        
        aip1ii = coords['a%d%d%d' % (nextJ,nextJ,j)]
        a0, a1 = aip1ii[0], aip1ii[1]
        qaip1ip1i = computeCircleForRegularHex(p, a0, a1)
        cip1ip1i = qaip1ip1i / func3(a0, a1)
        wip1ip1i = lambda x, y:  func33(x, y) * cip1ip1i #/ computeCircleForRegularHex(p, x, y) # * cip1ip1i
        # wip1ip1i = lambda x, y:  func4(x, y) * cip1ip1i

        return (wi,wiiip1, wip1ip1i,coeffs_2,coeffs_iiip1,coeffs_ip1ip1i)
    
    def computeFourthOrderFunction(self, j, p):
        coords,lines = self.getCoordsAndLines()
        line = lambda i: lines['d%d' %   hexCycle(i)]
        nextJ =   hexCycle(j+1)

        '''First terms of basis: wi'''

        #coefficients for the first function
        a_3 = sympy.Symbol('a_3_%d'%j)
        b_3 = sympy.Symbol('b_3_%d'%j)
        c_3 = sympy.Symbol('c_3_%d'%j)
        d_3 = sympy.Symbol('d_3_%d'%j)
        e_3 = sympy.Symbol('e_3_%d'%j)
        f_3 = sympy.Symbol('f_3_%d'%j)
        g_3 = sympy.Symbol('g_3_%d'%j)
        h_3 = sympy.Symbol('h_3_%d'%j)
        i_3 = sympy.Symbol('i_3_%d'%j)
        j_3 = sympy.Symbol('j_3_%d'%j)
        coeffs_3 = [a_3,b_3,c_3,d_3,e_3,f_3,g_3,h_3,i_3,j_3]

        func1 = lambda x,y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y)* computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) #* computeCubic(coeffs_3,x,y)
        func11 = lambda x, y: func1(x,y) * computeCubic(coeffs_3, x, y)
        a0, a1 = coords['a%d' % j][0], coords['a%d' % j][1]
        qai = computeCircleForRegularHex(p, a0, a1)
        ci = qai / func1(a0, a1)
        wi = lambda x, y: func11(x, y) * ci #/ computeCircleForRegularHex(p, x, y) # * ci
        # wi = lambda x, y:  func1(x, y) * ci

        '''Second terms of basis: wiiip1'''

        #coefficients for the second function
        a_iiip1 = sympy.Symbol('a_iiip1_%d'%j)
        b_iiip1 = sympy.Symbol('b_iiip1_%d'%j)
        c_iiip1 = sympy.Symbol('c_iiip1_%d'%j)
        d_iiip1 = sympy.Symbol('d_iiip1_%d'%j)
        e_iiip1 = sympy.Symbol('e_iiip1_%d'%j)
        f_iiip1 = sympy.Symbol('f_iiip1_%d'%j)

        coeffs_iiip1 = [a_iiip1, b_iiip1, c_iiip1, d_iiip1, e_iiip1, f_iiip1]
        func2 = lambda x,y:  computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) * computeLi(line(j), x, y) #* cls.computeConic(coeffs_iiip1, x,y)
        func22 = lambda x, y: func2(x,y) * computeConic(coeffs_iiip1, x, y)
        
        aiiip1 = coords['a%d%d%d' % (j,j,nextJ)]
        a0, a1 = aiiip1[0], aiiip1[1]
        qaiiip1 = computeCircleForRegularHex(p, a0, a1)
        ciiip1 = qaiiip1 / func2(a0, a1)
        wiiip1 = lambda x, y:  func22(x, y) * ciiip1 #/ computeCircleForRegularHex(p, x, y) # * ciiip1
        # wiiip1 = lambda x, y:  func2(x, y) * ciiip1


        '''Third terms of basis: wiip1'''

        #coefficients for the third function
   
        a_iip1 = sympy.Symbol('a_iip1_%d'%j)
        b_iip1 = sympy.Symbol('b_iip1_%d'%j)
        c_iip1 = sympy.Symbol('c_iip1_%d'%j)
        d_iip1 = sympy.Symbol('d_iip1_%d'%j)
        e_iip1 = sympy.Symbol('e_iip1_%d'%j)
        f_iip1 = sympy.Symbol('f_iip1_%d'%j)
        coeffs_iip1 = [a_iip1, b_iip1, c_iip1, d_iip1, e_iip1, f_iip1]
        
        func3 = lambda x,y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) * computeLi(line(j), x, y) #* computeCircleForRegularHex(p, x, y, 13./48)
        func33 = lambda x, y: func3(x,y) * computeConic(coeffs_iip1, x, y)
        ip1 = lambda i:   hexCycle(i)
        aiip1 = coords['a%d%d' % (j, ip1(j+1))]
        qaiip1 = computeCircleForRegularHex(p, aiip1[0], aiip1[1])
        ciip1 = qaiip1 / func3(aiip1[0], aiip1[1])
        wiip1 = lambda x, y:  func33(x, y) * ciip1 #/ computeCircleForRegularHex(p, x, y) # * ciip1
        # wiip1 = lambda x, y:  func3(x, y) * ciip1


        '''Fourth terms of basis: wip1ip1i'''

        #Coefficients for the fourth function
        a_ip1ip1i = sympy.Symbol('a_ip1ip1i_%d'%j)
        b_ip1ip1i = sympy.Symbol('b_ip1ip1i_%d'%j)
        c_ip1ip1i = sympy.Symbol('c_ip1ip1i_%d'%j)
        d_ip1ip1i = sympy.Symbol('d_ip1ip1i_%d'%j)
        e_ip1ip1i = sympy.Symbol('e_ip1ip1i_%d'%j)
        f_ip1ip1i = sympy.Symbol('f_ip1ip1i_%d'%j)

        coeffs_ip1ip1i = [a_ip1ip1i, b_ip1ip1i, c_ip1ip1i, d_ip1ip1i, e_ip1ip1i, f_ip1ip1i]
        func4 = lambda x,y:  computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) * computeLi(line(j), x, y) #* cls.computeConic(coeffs_ip1ip1i, x,y)
        func44 = lambda x, y: func4(x,y) * computeConic(coeffs_ip1ip1i, x, y)
        
        aip1ii = coords['a%d%d%d' % (nextJ,nextJ,j)]
        a0, a1 = aip1ii[0], aip1ii[1]
        qaip1ip1i = computeCircleForRegularHex(p, a0, a1)
        cip1ip1i = qaip1ip1i / func4(a0, a1)
        wip1ip1i = lambda x, y:  func44(x, y) * cip1ip1i #/ computeCircleForRegularHex(p, x, y) # * cip1ip1i
        # wip1ip1i = lambda x, y:  func4(x, y) * cip1ip1i

        return (wi,wiiip1, wiip1, wip1ip1i,coeffs_3,coeffs_iiip1, coeffs_iip1, coeffs_ip1ip1i)
    
    def computeFifthOrderFunction(self, j,p):
        coords,lines = self.getCoordsAndLines()
        line = lambda i: lines['d%d' %   hexCycle(i)]
        nextJ =   hexCycle(j+1)

        '''First terms of basis: wi'''

        #coefficients for the first function
        a_4 = sympy.Symbol('a_4_%d'%j)
        b_4 = sympy.Symbol('b_4_%d'%j)
        c_4 = sympy.Symbol('c_4_%d'%j)
        d_4 = sympy.Symbol('d_4_%d'%j)
        e_4 = sympy.Symbol('e_4_%d'%j)
        f_4 = sympy.Symbol('f_4_%d'%j)
        g_4 = sympy.Symbol('g_4_%d'%j)
        h_4 = sympy.Symbol('h_4_%d'%j)
        i_4 = sympy.Symbol('i_4_%d'%j)
        j_4 = sympy.Symbol('j_4_%d'%j)
        k_4 = sympy.Symbol('k_4_%d'%j)
        l_4 = sympy.Symbol('l_4_%d'%j)
        m_4 = sympy.Symbol('m_4_%d'%j)
        n_4 = sympy.Symbol('n_4_%d'%j)
        o_4 = sympy.Symbol('o_4_%d'%j)
        
        coeffs_4 = [a_4,b_4,c_4,d_4,e_4,f_4,g_4,h_4,i_4,j_4,k_4,l_4,m_4,n_4,o_4]

        func1 = lambda x,y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y)* computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) #* computeCubic(coeffs_3,x,y)
        func11 = lambda x, y: func1(x,y) * computeQuartic(coeffs_4, x, y)
        a0, a1 = coords['a%d' % j][0], coords['a%d' % j][1]
        qai = computeCircleForRegularHex(p, a0, a1)
        ci = qai / func1(a0, a1)
        wi = lambda x, y: func11(x, y) * ci #/ computeCircleForRegularHex(p, x, y) # * ci
        # wi = lambda x, y:  func1(x, y) * ci

        '''Second terms of basis: wi_1'''

        #coefficients for the second function
        a_3_1 = sympy.Symbol('a_3_1_%d'%j)
        b_3_1 = sympy.Symbol('b_3_1_%d'%j)
        c_3_1 = sympy.Symbol('c_3_1_%d'%j)
        d_3_1 = sympy.Symbol('d_3_1_%d'%j)
        e_3_1 = sympy.Symbol('e_3_1_%d'%j)
        f_3_1 = sympy.Symbol('f_3_1_%d'%j)
        g_3_1 = sympy.Symbol('g_3_1_%d'%j)
        h_3_1 = sympy.Symbol('h_3_1_%d'%j)
        i_3_1 = sympy.Symbol('i_3_1_%d'%j)
        j_3_1 = sympy.Symbol('j_3_1_%d'%j)
        
        coeffs_3_1 = [a_3_1, b_3_1, c_3_1, d_3_1, e_3_1, f_3_1, g_3_1, h_3_1, i_3_1, j_3_1]
        funcj_1 = lambda x,y:  computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) * computeLi(line(j), x, y) #* cls.computeConic(coeffs_iiip1, x,y)
        funcj_11 = lambda x, y: funcj_1(x,y) * computeCubic(coeffs_3_1, x, y)
        
        ai_j_1 = coords['a%d1' % j]
        a0, a1 = ai_j_1[0], ai_j_1[1]
        qai_j_1 = computeCircleForRegularHex(p, a0, a1)
        cij_1 = qai_j_1  / funcj_1(a0, a1)
        wi_1 = lambda x, y:  funcj_11(x, y) * cij_1 #/ computeCircleForRegularHex(p, x, y) # * ciiip1
        # wiiip1 = lambda x, y:  func2(x, y) * ciiip1


        '''Third terms of basis: wiip1'''

        #coefficients for the third function
        a_3_2 = sympy.Symbol('a_3_2_%d'%j)
        b_3_2 = sympy.Symbol('b_3_2_%d'%j)
        c_3_2 = sympy.Symbol('c_3_2_%d'%j)
        d_3_2 = sympy.Symbol('d_3_2_%d'%j)
        e_3_2 = sympy.Symbol('e_3_2_%d'%j)
        f_3_2 = sympy.Symbol('f_3_2_%d'%j)
        g_3_2 = sympy.Symbol('g_3_2_%d'%j)
        h_3_2 = sympy.Symbol('h_3_2_%d'%j)
        i_3_2 = sympy.Symbol('i_3_2_%d'%j)
        j_3_2 = sympy.Symbol('j_3_2_%d'%j)
        
        
        coeffs_3_2 = [a_3_2, b_3_2, c_3_2, d_3_2, e_3_2, f_3_2, g_3_2, h_3_2, i_3_2, j_3_2]
        funcj_2 = lambda x,y:  computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) * computeLi(line(j), x, y) #* cls.computeConic(coeffs_iiip1, x,y)
        funcj_22 = lambda x, y: funcj_2(x,y) * computeCubic(coeffs_3_2, x, y)
        
        ai_j_2 = coords['a%d2' % j]
        a0, a1 = ai_j_2[0], ai_j_2[1]
        qai_j_2 = computeCircleForRegularHex(p, a0, a1)
        cij_2 = qai_j_2 / funcj_2(a0, a1)
        wi_2 = lambda x, y:  funcj_22(x, y) * cij_2 

        '''Fourth terms of basis: wip1ip1i'''

        #Coefficients for the fourth function
        a_3_3 = sympy.Symbol('a_3_3_%d'%j)
        b_3_3 = sympy.Symbol('b_3_3_%d'%j)
        c_3_3 = sympy.Symbol('c_3_3_%d'%j)
        d_3_3 = sympy.Symbol('d_3_3_%d'%j)
        e_3_3 = sympy.Symbol('e_3_3_%d'%j)
        f_3_3 = sympy.Symbol('f_3_3_%d'%j)
        g_3_3 = sympy.Symbol('g_3_3_%d'%j)
        h_3_3 = sympy.Symbol('h_3_3_%d'%j)
        i_3_3 = sympy.Symbol('i_3_3_%d'%j)
        j_3_3 = sympy.Symbol('j_3_3_%d'%j)
        
        
        coeffs_3_3 = [a_3_3, b_3_3, c_3_3, d_3_3, e_3_3, f_3_3, g_3_3, h_3_3, i_3_3, j_3_3]
        funcj_3 = lambda x,y:  computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) * computeLi(line(j), x, y) #* cls.computeConic(coeffs_iiip1, x,y)
        funcj_33 = lambda x, y: funcj_3(x,y) * computeCubic(coeffs_3_3, x, y)
        
        ai_j_3 = coords['a%d3' % j]
        a0, a1 = ai_j_3[0], ai_j_3[1]
        qai_j_3 = computeCircleForRegularHex(p, a0, a1)
        cij_3 = qai_j_3 / funcj_3(a0, a1)
        wi_3 = lambda x, y:  funcj_33(x, y) * cij_3 
        
        '''Fifth terms of basis: wip1ip1i'''

        #Coefficients for the fifth function
        a_3_4 = sympy.Symbol('a_3_4_%d'%j)
        b_3_4 = sympy.Symbol('b_3_4_%d'%j)
        c_3_4 = sympy.Symbol('c_3_4_%d'%j)
        d_3_4 = sympy.Symbol('d_3_4_%d'%j)
        e_3_4 = sympy.Symbol('e_3_4_%d'%j)
        f_3_4 = sympy.Symbol('f_3_4_%d'%j)
        g_3_4 = sympy.Symbol('g_3_4_%d'%j)
        h_3_4 = sympy.Symbol('h_3_4_%d'%j)
        i_3_4 = sympy.Symbol('i_3_4_%d'%j)
        j_3_4 = sympy.Symbol('j_3_4_%d'%j)
        
        
        coeffs_3_4 = [a_3_4, b_3_4, c_3_4, d_3_4, e_3_4, f_3_4, g_3_4, h_3_4, i_3_4, j_3_4]
        funcj_4 = lambda x,y:  computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) * computeLi(line(j), x, y) #* cls.computeConic(coeffs_iiip1, x,y)
        funcj_44 = lambda x, y: funcj_4(x,y) * computeCubic(coeffs_3_4, x, y)
        
        ai_j_4 = coords['a%d4' % j]
        a0, a1 = ai_j_4[0], ai_j_4[1]
        qai_j_4 = computeCircleForRegularHex(p, a0, a1)
        cij_4 = qai_j_4 / funcj_4(a0, a1)
        wi_4 = lambda x, y:  funcj_44(x, y) * cij_4


        return (wi,wi_1, wi_2, wi_3,wi_4, coeffs_4, coeffs_3_1, coeffs_3_2, coeffs_3_3, coeffs_3_4)
    
    @classmethod
    def get_coeffs_Gout(cls, coords, lines,p):

        coeffs_2 = []
        coeffs_iiip1 = []
        coeffs_ip1ip1i = []

        line = lambda i: lines['d%d' %   hexCycle(i)]
        linep = lambda i, k: lines['dp%d%d' % (  hexCycle(i),   hexCycle(k))]      
        
        for j in range(1,7):
            coeffs_2_i = [1., 1.,0,0,0, -7./27*p**2]
            coeffs_iiip1_i = linep(j+1, j+3) #avant linep(j,j+4)
            coeffs_ip1ip1i_i = linep(j, j+4) #avant linep(j+1,j+3)
            nextJ = j + 1
    
            
            a0, a1 = coords['a%d'% j]
            if(j == 6):
                aiiip1 = coords['a%d%d%d' % (j,j,1)]
                aip1ip1i = coords['a%d%d%d' % (1,1,j)]
            else:
                aiiip1 = coords['a%d%d%d' % (j,j,nextJ)]
                aip1ip1i = coords['a%d%d%d' % (nextJ,nextJ,j)]
                
            aiiip1x, aiiip1y = aiiip1[0], aiiip1[1]
            aip1ip1ix, aip1ip1iy = aip1ip1i[0], aip1ip1i[1]
            
            func3 = lambda x,y:  computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) * computeLi(line(j), x, y) 
            
            ci = computeCircleForRegularHex(coeffs_2_i, a0, a1, False)
            ciiip1 = computeCircleForRegularHex(p, aiiip1x, aiiip1y)/(computeLi(linep(j+1, j+3), aiiip1x, aiiip1y)*func3(aiiip1x,aiiip1y))
            
            cip1ip1i = computeCircleForRegularHex(p, aip1ip1ix, aip1ip1iy)/ (computeLi(linep(j, j+4), aip1ip1ix, aip1ip1iy)* func3(aip1ip1ix, aip1ip1iy))
            
            for i in range(len(coeffs_2_i)):
                coeffs_2.append(coeffs_2_i[i]/ci)
            for h in range(len(coeffs_iiip1_i)):
                coeffs_iiip1.append(coeffs_iiip1_i[h]/ciiip1)
                coeffs_ip1ip1i.append(coeffs_ip1ip1i_i[h]/cip1ip1i)

        coeffs_Gout = coeffs_2 + coeffs_iiip1 + coeffs_ip1ip1i
        print("----------------------------------------------")
        print(func3(aiiip1x,aiiip1y))
        print(computeCircleForRegularHex(p, aiiip1x, aiiip1y))
        print(func3(aiiip1x,aiiip1y)/computeCircleForRegularHex(p, aiiip1x, aiiip1y))
        print(computeLi(coeffs_iiip1_i, aiiip1x,aiiip1y))
        print("----------------------------------------------")
        return coeffs_Gout
        
    


def plotFunctionOnHexagon(func, p, coords):
    maxX = p/(3**0.5)
    maxY = 0.5*p
    size = 750
    sizeX = 2*maxX/size
    sizeY = 2*maxY/size

    x, y = np.meshgrid(np.arange(-maxX, maxX, sizeX), np.arange(-maxY, maxY, sizeY))
    xf, yf = x.flatten(), y.flatten()

    points = np.vstack((xf, yf)).T
    path = Path(coords)
    flags = path.contains_points(points).reshape(size, size)
    it = np.nditer(flags, flags=['multi_index'])
    f = func(x, y)
    while not it.finished:
        if not flags[it.multi_index]:
            f[it.multi_index] = np.nan
        it.iternext()

    fig = plt.figure()
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    surf = ax.plot_surface(x, y, f)
    ax = fig.add_subplot(1, 2, 2)
    plt.contourf(x, y, f, 30)
    fig.colorbar(surf)
    plt.tight_layout()
    plt.show()
