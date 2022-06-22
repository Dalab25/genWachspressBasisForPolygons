# -*- coding: utf-8 -*-
#
# File: hexagonal_functions_generalised_basis.py
#
"""
This class defines the Wachspress shape functions by the general
methodology from D. Labeurthre. The elementary function coefficients 
for building the basis functions are stored in an ASCII file provided 
by the user. For orders 1 and 2, they are the same functions defined by Gout
in his PhD.
Reference: J. L. Gout, Elements finis polygonaux de Wachspress (1980)
"""
import numpy as np

from geometrical_elements import computeLi,computeCubic, computeConic, computeAdjointForIrregularPentagon, pentaCycle
from WSF_to_generate_irregular_pentagon import WSFTogenerateIrregularPentagon

class PentagonalFunctionsGeneralisedBasis(object):

    def __init__(self, order, pentagon, p, coeffsFile, homothetic=False):
        """
        Initialises the calculator for computing the functions at a
        given order using user-defined coefficients.
        """
        self.__pentagon = pentagon
        self.__order = order
        self.__pitch = p
        self.__coords, self.__lines = WSFTogenerateIrregularPentagon(self.__order, self.__pentagon).coords, WSFTogenerateIrregularPentagon(self.__order, self.__pentagon).lines
        
        if self.__order > 2:
            self.__coeffs = np.loadtxt(coeffsFile)
            # as from order n > 2, for each side, this method requires
            # n(n+1)/2 coefficients to define the polynomial function
            # for the WSF anchored at the vertex of the hexagon. For
            # functions which are anchored at the remaining nodes on
            # the edge of the hexagon, it is 1 degree lower than, such
            # that the number of coefficients is (n-1)(1-1+1)/2 =
            # (n-1)n/2, and there are (n-1) such nodes. Thus, the
            # total number of coeffcients per side is
            # n(n+1)/2 + (n-1)*(n-1)/2
            assert(self.__coeffs.size == (self.__order*(self.__order+1)//2 + (self.__order-1) * (self.__order-1)*self.__order//2)*5)

        if homothetic:
            self.__homFact = 1./self.__pitch
        else:
            self.__homFact = 1.
            
    @property 
    def pentagon(self):
        return self.__pentagon
    
    @property
    def coords(self):
        return self.__coords

    @property
    def lines(self):
        return self.__lines

    def getCoordsLines(self):
        return self.__coords, self.__lines

    def computeFirstOrderFunction(self, j):
        """
        Function to build first-order Wachspress elements (cf. Gout 1985).
        """
        return WSFTogenerateIrregularPentagon(self.__order, self.__pentagon).computeFirstOrderFunction(j)

    def computeSecondOrderFunction(self, j):
        """
        Function to build second-order Wachspress elements (cf. Gout 1985).
        """
        return WSFTogenerateIrregularPentagon(self.__order, self.__pentagon).computeSecondOrderFunction(j)

    def computeThirdOrderFunction(self, j):
        """
        Function to build third-order Wachspress elements with user-defined coefficients.
        """
        # First terms of basis: wi
        line = lambda i: self.__lines['d%d' % pentaCycle(i)]
        nextJ = pentaCycle(j+1)

        # s represents the number of coefficient for the function ri
        s = self.__order*(self.__order+1)//2

        coeffs_2 = self.__coeffs[s*j:s*(j+1)]
        func1 = lambda x, y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) 
        func11 = lambda x, y: computeConic(coeffs_2, x, y)

        a0, a1 = self.__coords['a%d' % j][0], self.__coords['a%d' % j][1]
        qai = computeAdjointForIrregularPentagon(self.__pentagon.coeffsAdjoint, a0*self.__homFact, a1*self.__homFact)
        ci = qai / func1(a0, a1)

        wi = lambda x, y: func1(x, y) * func11(x*self.__homFact, y*self.__homFact) * ci / computeAdjointForIrregularPentagon(self.__pentagon.coeffsAdjoint, x*self.__homFact, y*self.__homFact)

        # Second terms of basis: wiiip1 coefficients for the second function
        # s_1 is the number of coefficient for the functions rij
        # we need to multiply s by 5 because there are 5 functions wi
        s *= 5
        s_1 = self.__order*(self.__order-1)//2

        coeffs_iiip1 = self.__coeffs[s + s_1*(j):s + s_1*(j+1)]
        func2 = lambda x, y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j), x, y)
        func22 = lambda x, y: computeLi(coeffs_iiip1, x, y)

        aiiip1 = self.__coords['a%d%d%d' % (j, j, nextJ)]
        a0, a1 = aiiip1[0], aiiip1[1]
        qaiiip1 = computeAdjointForIrregularPentagon(self.__pentagon.coeffsAdjoint, a0*self.__homFact, a1*self.__homFact)
        ciiip1 = qaiiip1 / func2(a0, a1)

        wiiip1 = lambda x, y:  func2(x, y) * func22(x*self.__homFact, y*self.__homFact) * ciiip1 / computeAdjointForIrregularPentagon(self.__pentagon.coeffsAdjoint, x*self.__homFact, y*self.__homFact)

        # Third terms of basis: wip1ip1i coefficients for the third function
        # s_2 is the number of coefficient for the functions rij
        # we need to multiply s_1 by 5 because there is 6 functions wi_1
        # to use the right coefficients corresponding to wi_2 we need
        # to do a sum, see line s = s + 5*s_1
        s_2 = s_1
        s = s + 5*s_1

        coeffs_ip1ip1i = self.__coeffs[s + s_2*j:s + s_2*(j+1)]

        func3 = lambda x, y:  computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j), x, y)
        func33 = lambda x, y: computeLi(coeffs_ip1ip1i, x, y)

        aip1ii = self.__coords['a%d%d%d' % (nextJ, nextJ, j)]
        a0, a1 = aip1ii[0], aip1ii[1]
        qaip1ip1i = computeAdjointForIrregularPentagon(self.__pentagon.coeffsAdjoint, a0*self.__homFact, a1*self.__homFact)
        cip1ip1i = qaip1ip1i / func3(a0, a1)

        wip1ip1i = lambda x, y: func3(x, y) * func33(x*self.__homFact, y*self.__homFact) * cip1ip1i/ computeAdjointForIrregularPentagon(self.__pentagon.coeffsAdjoint, x*self.__homFact, y*self.__homFact)

        return (wi, wiiip1, wip1ip1i)

    def computeFourthOrderFunction(self, j):
        """
        Function to build fourth-order Wachspress elements with
        user-defined coefficients.
        """
        # First terms of basis: wi
        line = lambda i: self.__lines['d%d' % pentaCycle(i)]
        nextJ = pentaCycle(j+1)

        s = self.__order*(self.__order+1)//2
        coeffs_3 = self.__coeffs[s*j:s*(j+1)]
        func1 = lambda x, y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) 
        func11 = lambda x, y: computeCubic(coeffs_3, x, y)

        a0, a1 = self.__coords['a%d' % j][0], self.__coords['a%d' % j][1]
        qai = computeAdjointForIrregularPentagon(self.__pentagon.coeffsAdjoint, a0*self.__homFact, a1*self.__homFact)
        ci = qai / func1(a0, a1)

        wi = lambda x, y: func1(x, y) * func11(x*self.__homFact, y*self.__homFact) * ci / computeAdjointForIrregularPentagon(self.__pentagon.coeffsAdjoint, x*self.__homFact, y*self.__homFact)

        # Second terms of basis: wiiip1
        s *= 5
        s_1 = self.__order*(self.__order-1)//2

        # coefficients for the second function
        coeffs_iiip1 = self.__coeffs[s + s_1*(j):s + s_1*(j+1)]
        func2 = lambda x, y:  computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j), x, y)
        func22 = lambda x, y: computeConic(coeffs_iiip1, x, y)

        aiiip1 = self.__coords['a%d%d%d' % (j, j, nextJ)]
        a0, a1 = aiiip1[0], aiiip1[1]
        qaiiip1 = computeAdjointForIrregularPentagon(self.__pentagon.coeffsAdjoint, a0*self.__homFact, a1*self.__homFact)
        ciiip1 = qaiiip1 / func2(a0, a1)

        wiiip1 = lambda x, y:  func2(x, y) * func22(x*self.__homFact, y*self.__homFact) * ciiip1 / computeAdjointForIrregularPentagon(self.__pentagon.coeffsAdjoint, x*self.__homFact, y*self.__homFact)

        # Third terms of basis: wiip1
        s_2 = s_1
        s = s + 5*s_1

        coeffs_iip1 = self.__coeffs[s + s_2*(j):s + s_2*(j+1)]

        func3 = lambda x, y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) *  computeLi(line(j), x, y)
        func33 = lambda x, y: computeConic(coeffs_iip1, x, y)

        ip1 = lambda i: pentaCycle(i)
        aiip1 = self.__coords['a%d%d' % (j, ip1(j+1))]
        qaiip1 = computeAdjointForIrregularPentagon(self.__pentagon.coeffsAdjoint, aiip1[0]*self.__homFact, aiip1[1]*self.__homFact)
        ciip1 = qaiip1 / func3(aiip1[0], aiip1[1])

        wiip1 = lambda x, y:  func3(x, y) * func33(x*self.__homFact, y*self.__homFact) * ciip1 / computeAdjointForIrregularPentagon(self.__pentagon.coeffsAdjoint,  x*self.__homFact, y*self.__homFact)

        # Fourth terms of basis: wip1ip1i
        s_3 = s_1
        s = s + 5*s_1

        # Coefficients for the fourth function
        coeffs_ip1ip1i = self.__coeffs[s + s_3*(j):s + s_3*(j+1)]

        func4 = lambda x, y:  computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j), x, y)
        func44 = lambda x, y: computeConic(coeffs_ip1ip1i, x, y)

        aip1ii = self.__coords['a%d%d%d' % (nextJ, nextJ, j)]
        a0, a1 = aip1ii[0], aip1ii[1]
        qaip1ip1i = computeAdjointForIrregularPentagon(self.__pentagon.coeffsAdjoint, a0*self.__homFact, a1*self.__homFact)
        cip1ip1i = qaip1ip1i / func4(a0, a1)

        wip1ip1i = lambda x, y:  func4(x, y) * func44(x*self.__homFact, y*self.__homFact) * cip1ip1i  / computeAdjointForIrregularPentagon(self.__pentagon.coeffsAdjoint, x*self.__homFact, y*self.__homFact)

        return (wi, wiiip1, wiip1, wip1ip1i)
    