# -*- coding: utf-8 -*-
#
# File: pentagonal_functions_generalised_basis.py
#
"""
This class defines the Wachspress shape functions by the general
methodology from D. Labeurthre.
The elementary function coefficients for building the basis functions
are stored in an ASCII file provided by the user.
For orders 1 and 2, they are the same functions defined by Gout.
Reference: J. L. Gout, Rational Wachspress-type Finite Elements on
Regular Hexagons, Journal of Numerical Analysis (1985)
"""
import numpy as np

from geometrical_elements import computeLi, computeQuartic, computeCubic, computeConic, computeCircleForRegularPolygon, hexCycle
from hexagonal_nodal_dofs import buildCoordsAndLinesForWSF
from hexagonal_functions_gout_basis import HexagonalFunctionsGoutBasis


class HexagonalFunctionsGeneralisedBasis(object):

    def __init__(self, hexMesh, order, coeffsFile=None, homothetic=True):
        """
        Initialises the calculator for computing the functions at a
        given order using user-defined coefficients.
        """
        self.__hexMesh = hexMesh
        self.__order = order
        self.__pitch = self.__hexMesh.pitch
        self.__coords, self.__lines = buildCoordsAndLinesForWSF(hexMesh.halfPitch,
                                                                    hexMesh.halfSide,
                                                                    self.__order
                                                                )

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
            assert(self.__coeffs.size == (self.__order*(self.__order+1)//2 + (self.__order-1) * (self.__order-1)*self.__order//2)*6)

        if homothetic:
            self.__homFact = 1./self.__pitch
        else:
            self.__homFact = 1.

    def getCoordsLines(self):
        return self.__coords, self.__lines

    def computeFirstOrderFunction(self, j):
        """
        Function to build first-order Wachspress elements (cf. Gout 1985).
        """
        return HexagonalFunctionsGoutBasis(self.__hexMesh, self.__order).computeFirstOrderFunction(j)

    def computeSecondOrderFunction(self, j):
        """
        Function to build second-order Wachspress elements (cf. Gout 1985).
        """
        return HexagonalFunctionsGoutBasis(self.__hexMesh, self.__order).computeSecondOrderFunction(j)

    def computeThirdOrderFunction(self, j):
        """
        Function to build third-order Wachspress elements with user-defined coefficients.
        """
        # First terms of basis: wi
        line = lambda i: self.__lines['d%d' % hexCycle(i)]
        nextJ = hexCycle(j+1)

        # s represents the number of coefficient for the function ri
        s = self.__order*(self.__order+1)//2

        coeffs_2 = self.__coeffs[s*(j - 1):s*j]
        func1 = lambda x, y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y)
        func11 = lambda x, y: computeConic(coeffs_2, x, y)

        a0, a1 = self.__coords['a%d' % j][0], self.__coords['a%d' % j][1]
        qai = computeCircleForRegularPolygon(self.__pitch, a0*self.__homFact, a1*self.__homFact, self.__homFact**2)
        ci = qai / func1(a0, a1)

        wi = lambda x, y: func1(x, y) * func11(x*self.__homFact, y*self.__homFact) * ci / computeCircleForRegularPolygon(self.__pitch, x*self.__homFact, y*self.__homFact, self.__homFact**2)

        # Second terms of basis: wiiip1 coefficients for the second function
        # s_1 is the number of coefficient for the functions rij
        # we need to multiply s by 6 because there are 6 functions wi
        s *= 6
        s_1 = self.__order*(self.__order-1)//2

        coeffs_iiip1 = self.__coeffs[s + s_1*(j-1):s + s_1*j]
        func2 = lambda x, y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) * computeLi(line(j), x, y)
        func22 = lambda x, y: computeLi(coeffs_iiip1, x, y)

        aiiip1 = self.__coords['a%d%d%d' % (j, j, nextJ)]
        a0, a1 = aiiip1[0], aiiip1[1]
        qaiiip1 = computeCircleForRegularPolygon(self.__pitch, a0*self.__homFact, a1*self.__homFact, self.__homFact**2)
        ciiip1 = qaiiip1 / func2(a0, a1)

        wiiip1 = lambda x, y:  func2(x, y) * func22(x*self.__homFact, y*self.__homFact) * ciiip1 / computeCircleForRegularPolygon(self.__pitch, x*self.__homFact, y*self.__homFact, self.__homFact**2)

        # Third terms of basis: wip1ip1i coefficients for the third function
        # s_2 is the number of coefficient for the functions rij
        # we need to multiply s_1 by 6 because there is 6 functions wi_1
        # to use the right coefficients corresponding to wi_2 we need
        # to do a sum, see line s = s + 6*s_1
        s_2 = s_1
        s = s + 6*s_1

        coeffs_ip1ip1i = self.__coeffs[s + s_2*(j-1):s + s_2*j]

        func3 = lambda x, y:  computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) * computeLi(line(j), x, y)
        func33 = lambda x, y: computeLi(coeffs_ip1ip1i, x, y)

        aip1ii = self.__coords['a%d%d%d' % (nextJ, nextJ, j)]
        a0, a1 = aip1ii[0], aip1ii[1]
        qaip1ip1i = computeCircleForRegularPolygon(self.__pitch, a0*self.__homFact, a1*self.__homFact, self.__homFact**2)
        cip1ip1i = qaip1ip1i / func3(a0, a1)

        wip1ip1i = lambda x, y: func3(x, y) * func33(x*self.__homFact, y*self.__homFact) * cip1ip1i / computeCircleForRegularPolygon(self.__pitch, x*self.__homFact, y*self.__homFact, self.__homFact**2)

        return (wi, wiiip1, wip1ip1i)

    def computeFourthOrderFunction(self, j):
        """
        Function to build fourth-order Wachspress elements with
        user-defined coefficients.
        """
        # First terms of basis: wi
        line = lambda i: self.__lines['d%d' % hexCycle(i)]
        nextJ = hexCycle(j+1)

        s = self.__order*(self.__order+1)//2
        coeffs_3 = self.__coeffs[s*(j - 1):s*j]
        func1 = lambda x, y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y)
        func11 = lambda x, y: computeCubic(coeffs_3, x, y)

        a0, a1 = self.__coords['a%d' % j][0], self.__coords['a%d' % j][1]
        qai = computeCircleForRegularPolygon(self.__pitch, a0*self.__homFact, a1*self.__homFact, self.__homFact**2)
        ci = qai / func1(a0, a1)

        wi = lambda x, y: func1(x, y) * func11(x*self.__homFact, y*self.__homFact) * ci / computeCircleForRegularPolygon(self.__pitch, x*self.__homFact, y*self.__homFact, self.__homFact**2)

        # Second terms of basis: wiiip1
        s *= 6
        s_1 = self.__order*(self.__order-1)//2

        # coefficients for the second function
        coeffs_iiip1 = self.__coeffs[s + s_1*(j-1):s + s_1*j]
        func2 = lambda x, y:  computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) * computeLi(line(j), x, y)
        func22 = lambda x, y: computeConic(coeffs_iiip1, x, y)

        aiiip1 = self.__coords['a%d%d%d' % (j, j, nextJ)]
        a0, a1 = aiiip1[0], aiiip1[1]
        qaiiip1 = computeCircleForRegularPolygon(self.__pitch, a0*self.__homFact, a1*self.__homFact, self.__homFact**2)
        ciiip1 = qaiiip1 / func2(a0, a1)

        wiiip1 = lambda x, y:  func2(x, y) * func22(x*self.__homFact, y*self.__homFact) * ciiip1 / computeCircleForRegularPolygon(self.__pitch, x*self.__homFact, y*self.__homFact, self.__homFact**2)

        # Third terms of basis: wiip1
        s_2 = s_1
        s = s + 6*s_1

        coeffs_iip1 = self.__coeffs[s + s_2*(j-1):s + s_2*j]

        func3 = lambda x, y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) * computeLi(line(j), x, y)
        func33 = lambda x, y: computeConic(coeffs_iip1, x, y)

        ip1 = lambda i: hexCycle(i)
        aiip1 = self.__coords['a%d%d' % (j, ip1(j+1))]
        qaiip1 = computeCircleForRegularPolygon(self.__pitch, aiip1[0]*self.__homFact, aiip1[1]*self.__homFact, self.__homFact**2)
        ciip1 = qaiip1 / func3(aiip1[0], aiip1[1])

        wiip1 = lambda x, y:  func3(x, y) * func33(x*self.__homFact, y*self.__homFact) * ciip1 / computeCircleForRegularPolygon(self.__pitch, x*self.__homFact, y*self.__homFact, self.__homFact**2)

        # Fourth terms of basis: wip1ip1i
        s_3 = s_1
        s = s + 6*s_1

        # Coefficients for the fourth function
        coeffs_ip1ip1i = self.__coeffs[s + s_3*(j-1):s + s_3*j]

        func4 = lambda x, y:  computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) * computeLi(line(j), x, y)
        func44 = lambda x, y: computeConic(coeffs_ip1ip1i, x, y)

        aip1ii = self.__coords['a%d%d%d' % (nextJ, nextJ, j)]
        a0, a1 = aip1ii[0], aip1ii[1]
        qaip1ip1i = computeCircleForRegularPolygon(self.__pitch, a0*self.__homFact, a1*self.__homFact, self.__homFact**2)
        cip1ip1i = qaip1ip1i / func4(a0, a1)

        wip1ip1i = lambda x, y:  func4(x, y) * func44(x*self.__homFact, y*self.__homFact) * cip1ip1i / computeCircleForRegularPolygon(self.__pitch, x*self.__homFact, y*self.__homFact, self.__homFact**2)

        return (wi, wiiip1, wiip1, wip1ip1i)

    def computeFifthOrderFunction(self, j):
        """
        Function to build fifth-order Wachspress elements with
        user-defined coefficients.
        """
        # First terms of basis: wi
        line = lambda i: self.__lines['d%d' % hexCycle(i)]
        nextJ = hexCycle(j+1)

        s = self.__order*(self.__order+1)//2
        # coefficients for the first function
        coeffs_4 = self.__coeffs[s*(j - 1):s*j]
        func1 = lambda x, y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y)
        func11 = lambda x, y: computeQuartic(coeffs_4, x, y)
        a0, a1 = self.__coords['a%d' % j][0], self.__coords['a%d' % j][1]
        qai = computeCircleForRegularPolygon(self.__pitch, a0*self.__homFact, a1*self.__homFact, self.__homFact**2)

        ci = qai / func1(a0, a1)

        wi = lambda x, y: func1(x, y)*func11(x*self.__homFact, y*self.__homFact) * ci / computeCircleForRegularPolygon(self.__pitch, x*self.__homFact, y*self.__homFact, self.__homFact**2)

        # Second terms of basis: wi_1
        s *= 6
        s_1 = self.__order*(self.__order-1)//2

        # coefficients for the second function
        coeffs_3_1 = self.__coeffs[s + s_1*(j-1):s + s_1*j]
        funcj_1 = lambda x, y:  computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) * computeLi(line(j), x, y)
        funcj_11 = lambda x, y: computeCubic(coeffs_3_1, x, y)

        ai_j_1 = self.__coords['a%d1' % j]
        a0, a1 = ai_j_1[0], ai_j_1[1]
        qai_j_1 = computeCircleForRegularPolygon(self.__pitch, a0*self.__homFact, a1*self.__homFact, self.__homFact**2)
        cij_1 = qai_j_1 / funcj_1(a0, a1)

        wi_1 = lambda x, y:  funcj_1(x, y) * funcj_11(x*self.__homFact, y*self.__homFact) * cij_1 / computeCircleForRegularPolygon(self.__pitch, x*self.__homFact, y*self.__homFact, self.__homFact**2)

        # Third terms of basis: wiip1
        s_2 = s_1
        s = s + 6*s_1

        # coefficients for the third function
        coeffs_3_2 = self.__coeffs[s + s_2*(j-1):s + s_2*j]

        funcj_2 = lambda x, y:  computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) * computeLi(line(j), x, y)
        funcj_22 = lambda x, y: computeCubic(coeffs_3_2, x, y)

        ai_j_2 = self.__coords['a%d2' % j]
        a0, a1 = ai_j_2[0], ai_j_2[1]
        qai_j_2 = computeCircleForRegularPolygon(self.__pitch, a0*self.__homFact, a1*self.__homFact, self.__homFact**2)
        cij_2 = qai_j_2 / funcj_2(a0, a1)

        wi_2 = lambda x, y:  funcj_2(x, y) * funcj_22(x*self.__homFact, y*self.__homFact) * cij_2 / computeCircleForRegularPolygon(self.__pitch, x*self.__homFact, y*self.__homFact, self.__homFact**2)

        # Fourth terms of basis: wip1ip1i
        s_3 = s_1
        s = s + 6*s_1

        coeffs_3_3 = self.__coeffs[s + s_3*(j-1):s + s_3*j]
        funcj_3 = lambda x, y:  computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) * computeLi(line(j), x, y)
        funcj_33 = lambda x, y: computeCubic(coeffs_3_3, x, y)

        ai_j_3 = self.__coords['a%d3' % j]
        a0, a1 = ai_j_3[0], ai_j_3[1]

        qai_j_3 = computeCircleForRegularPolygon(self.__pitch, a0*self.__homFact, a1*self.__homFact, self.__homFact**2)
        cij_3 = qai_j_3 / funcj_3(a0, a1)

        wi_3 = lambda x, y:  funcj_3(x, y) * funcj_33(x*self.__homFact, y*self.__homFact) * cij_3 / computeCircleForRegularPolygon(self.__pitch, x*self.__homFact, y*self.__homFact, self.__homFact**2)

        # Fifth terms of basis: wip1ip1i
        s_4 = s_1
        s = s + 6*s_1

        # Coefficients for the fifth function
        coeffs_3_4 = self.__coeffs[s + s_4*(j-1):s + s_4*j]
        funcj_4 = lambda x, y:  computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) * computeLi(line(j), x, y)
        funcj_44 = lambda x, y:  computeCubic(coeffs_3_4, x, y)

        ai_j_4 = self.__coords['a%d4' % j]
        a0, a1 = ai_j_4[0], ai_j_4[1]
        qai_j_4 = computeCircleForRegularPolygon(self.__pitch, a0*self.__homFact, a1*self.__homFact, self.__homFact**2)
        cij_4 = qai_j_4 / funcj_4(a0, a1)

        wi_4 = lambda x, y:  funcj_4(x, y) * funcj_44(x*self.__homFact, y*self.__homFact) * cij_4 / computeCircleForRegularPolygon(self.__pitch, x*self.__homFact, y*self.__homFact, self.__homFact**2)

        return (wi, wi_1, wi_2, wi_3, wi_4)
