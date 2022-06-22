# -*- coding: utf-8 -*-
#
# File: wachspress_shape_functions.py
"""
This class defines the Wachspress shape functions and the geometric
elements required for their construction.
To compute with the symbolic calulcation the basis, the denominator
is not taking into account.
"""
import sympy
from pentagonal_nodal_dofs import buildCoordsAndLinesForWSF
from geometrical_elements import  pentaCycle, nCycle, computeConic, computeCubic, computeLi, computeAdjointForIrregularPentagon


class WSFTogenerateIrregularPentagon(object):

    def __init__(self, order, pentagon):
        """ Constructor with a given order and hexMesh, required to
        obtain geometrical data such as lattice pitch, halfside of
        hexagon. It is important to note that the functions at 
        orders 1-2 are unique.
        """
        self.__order = order
        self.__pentagon = pentagon
        self.__basisSize = 5 * self.__order
        self.__coordsDict, self.__lines = buildCoordsAndLinesForWSF(pentagon, order)

        
        if self.__order == 1:
            self.__func = self.computeFirstOrderFunction
        elif self.__order == 2:
            self.__func = self.computeSecondOrderFunction
        elif self.__order == 3:
            self.__func = self.computeThirdOrderFunction
        elif self.__order == 4:
            self.__func = self.computeFourthOrderFunction
        else:
            raise ValueError("Order %s is not yet available" % order)

    @property
    def order(self):
        return self.__order

    @property
    def basisSize(self):
        return self.__basisSize

    @property
    def pentagon(self):
        return self.__pentagon
    
    @property
    def coords(self):
        return self.__coordsDict

    @property
    def lines(self):
        return self.__lines

    def computeFirstOrderFunction(self, j):
        """
        Function to build first-order Wachspress elements on a pentagon
        """
        line = lambda i: self.__lines['d%d' % nCycle(i, len(self.__pentagon.vertices))]

        func = lambda x, y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y)
        ajx, ajy = self.coords['a%d' % j][0], self.coords['a%d' % j][1]
        qai = computeAdjointForIrregularPentagon(self.__pentagon.coeffsAdjoint,ajx, ajy)
        ci = qai / func(ajx, ajy)
        wi = lambda x, y: ci * func(x, y) / computeAdjointForIrregularPentagon(self.__pentagon.coeffsAdjoint, x, y)
        return (wi)
    
    def computeSecondOrderFunction(self,j):
        """
        Function to build second-order Wachspress elements on a pentagon
        """
        #line = lambda i: self.__lines['d%d' % nCycle(i, len(self.__pentagon.vertices))]
        line = lambda i: self.lines['d%d' % pentaCycle(i)]
        linep = self.lines['dp%d' % j]
        
        func1 = lambda x, y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(linep, x, y)
        ajx, ajy = self.coords['a%d' % j][0], self.coords['a%d' % j][1]
        qai = computeAdjointForIrregularPentagon(self.__pentagon.coeffsAdjoint, ajx, ajy)
        ci = qai / func1(ajx, ajy)
        wi = lambda x, y: ci * func1(x, y) / computeAdjointForIrregularPentagon(self.__pentagon.coeffsAdjoint,x, y)

        func2 = lambda x, y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j), x, y)
        ip1 = lambda i: pentaCycle(i)
        aiip1 = self.coords['a%d%d' % (j, ip1(j+1))]
        qaiip1 = computeAdjointForIrregularPentagon(self.__pentagon.coeffsAdjoint, aiip1[0], aiip1[1])
        ciip1 = qaiip1 / func2(aiip1[0], aiip1[1])
        wiip1 = lambda x, y: ciip1 * func2(x, y) / computeAdjointForIrregularPentagon(self.__pentagon.coeffsAdjoint, x, y)

        return (wi,wiip1)

    def computeThirdOrderFunction(self,j):
        """
        Function to build third-order Wachspress elements on a pentagon
        """ 
        line = lambda i: self.lines['d%d' % pentaCycle(i)]
        nextJ = pentaCycle(j+1)
        #coefficients for the first function
        a_2 = sympy.Symbol('a_2_%d'%j)
        b_2 = sympy.Symbol('b_2_%d'%j)
        c_2 = sympy.Symbol('c_2_%d'%j)
        d_2 = sympy.Symbol('d_2_%d'%j)
        e_2 = sympy.Symbol('e_2_%d'%j)
        f_2 = sympy.Symbol('f_2_%d'%j)

        coeffs_2 = [a_2,b_2,c_2,d_2,e_2,f_2]

        func1 = lambda x,y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y)* computeLi(line(j+4), x, y) 
        func11 = lambda x, y: func1(x,y) * computeConic(coeffs_2, x, y)
        ajx, ajy = self.coords['a%d' % j][0], self.coords['a%d' % j][1]
        qai = computeAdjointForIrregularPentagon(self.__pentagon.coeffsAdjoint, ajx, ajy)
        ci = qai / func1(ajx, ajy)
        wi = lambda x, y: func11(x, y) * ci #/ computeCircleForRegularHex(p, x, y) # * ci
        # wi = lambda x, y:  func1(x, y) * ci

        '''Second terms of basis: wiiip1'''

        #coefficients for the second function
        a_iiip1 = sympy.Symbol('a_iiip1_%d'%j)
        b_iiip1 = sympy.Symbol('b_iiip1_%d'%j)
        c_iiip1 = sympy.Symbol('c_iiip1_%d'%j)

        coeffs_iiip1 = [a_iiip1, b_iiip1, c_iiip1]
        func2 = lambda x,y:  computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j), x, y) 
        func22 = lambda x, y: func2(x,y) * computeLi(coeffs_iiip1, x, y)
        
        aiiip1 = self.coords['a%d%d%d' % (j,j,nextJ)]
        a0, a1 = aiiip1[0], aiiip1[1]
        qaiiip1 = computeAdjointForIrregularPentagon(self.__pentagon.coeffsAdjoint, a0, a1)
        ciiip1 = qaiiip1 / func2(a0, a1)
        wiiip1 = lambda x, y:  func22(x, y) * ciiip1 #/ computeCircleForRegularHex(p, x, y) # * ciiip1
        # wiiip1 = lambda x, y:  func2(x, y) * ciiip1

        '''Third terms of basis: wip1ip1i'''

        #Coefficients for the third function
        a_ip1ip1i = sympy.Symbol('a_ip1ip1i_%d'%j)
        b_ip1ip1i = sympy.Symbol('b_ip1ip1i_%d'%j)
        c_ip1ip1i = sympy.Symbol('c_ip1ip1i_%d'%j)


        coeffs_ip1ip1i = [a_ip1ip1i, b_ip1ip1i, c_ip1ip1i]
        func3 = lambda x,y:  computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j), x, y)
        func33 = lambda x, y: func3(x,y) * computeLi(coeffs_ip1ip1i, x, y)
        
        aip1ii = self.coords['a%d%d%d' % (nextJ,nextJ,j)]
        a0, a1 = aip1ii[0], aip1ii[1]
        qaip1ip1i = computeAdjointForIrregularPentagon(self.__pentagon.coeffsAdjoint, a0, a1)
        cip1ip1i = qaip1ip1i / func3(a0, a1)
        wip1ip1i = lambda x, y:  func33(x, y) * cip1ip1i #/ computeCircleForRegularHex(p, x, y) # * cip1ip1i
        # wip1ip1i = lambda x, y:  func4(x, y) * cip1ip1i

        return (wi,wiiip1, wip1ip1i,coeffs_2,coeffs_iiip1,coeffs_ip1ip1i)

    def computeFourthOrderFunction(self, j):
    
    
            line = lambda i: self.lines['d%d' % pentaCycle(i)]
            nextJ = pentaCycle(j+1)
    
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
    
            func1 = lambda x,y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y)* computeLi(line(j+4), x, y)#* computeCubic(coeffs_3,x,y)
            func11 = lambda x, y: func1(x,y) * computeCubic(coeffs_3, x, y)
            a0, a1 = self.coords['a%d' % j][0],self. coords['a%d' % j][1]
            qai = computeAdjointForIrregularPentagon(self.__pentagon.coeffsAdjoint, a0, a1)
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
            func2 = lambda x,y:  computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j), x, y) #* computeConic(coeffs_iiip1, x,y)
            func22 = lambda x, y: func2(x,y) * computeConic(coeffs_iiip1, x, y)
            
            aiiip1 = self.coords['a%d%d%d' % (j,j,nextJ)]
            a0, a1 = aiiip1[0], aiiip1[1]
            qaiiip1 = computeAdjointForIrregularPentagon(self.__pentagon.coeffsAdjoint, a0, a1)
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
            
            func3 = lambda x,y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) *  computeLi(line(j), x, y) #* computeCircleForRegularHex(p, x, y, 13./48)
            func33 = lambda x, y: func3(x,y) * computeConic(coeffs_iip1, x, y)
            ip1 = lambda i: pentaCycle(i)
            aiip1 = self.coords['a%d%d' % (j, ip1(j+1))]
            qaiip1 = computeAdjointForIrregularPentagon(self.__pentagon.coeffsAdjoint, aiip1[0], aiip1[1])
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
            func4 = lambda x,y:  computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j), x, y) #* computeConic(coeffs_ip1ip1i, x,y)
            func44 = lambda x, y: func4(x,y) * computeConic(coeffs_ip1ip1i, x, y)
            
            aip1ii = self.coords['a%d%d%d' % (nextJ,nextJ,j)]
            a0, a1 = aip1ii[0], aip1ii[1]
            qaip1ip1i = computeAdjointForIrregularPentagon(self.__pentagon.coeffsAdjoint, a0, a1)
            cip1ip1i = qaip1ip1i / func4(a0, a1)
            wip1ip1i = lambda x, y:  func44(x, y) * cip1ip1i #/ computeCircleForRegularHex(p, x, y) # * cip1ip1i
            # wip1ip1i = lambda x, y:  func4(x, y) * cip1ip1i
    
            return (wi,wiiip1, wiip1, wip1ip1i,coeffs_3,coeffs_iiip1, coeffs_iip1, coeffs_ip1ip1i)