# -*- coding: utf-8 -*-
#
# File: hexagonal_functions_gout.py
#
"""
This class defines the Wachspress shape functions as given by Gout.
Reference: J. L. Gout, Rational Wachspress-type Finite Elements on
Regular Hexagons, Journal of Numerical Analysis (1985)
"""
from geometrical_elements import computeLi, computeCircleForRegularPolygon, hexCycle
from hexagonal_nodal_dofs import buildCoordsAndLinesForWSF


class HexagonalFunctionsGoutBasis(object):

    def __init__(self, hexMesh, order):
        """
        Initialises the calculator for computing the functions at a
        given order using Gout's method.
        """
        self.__pitch = hexMesh.pitch
        self.__coords, self.__lines = buildCoordsAndLinesForWSF(hexMesh.halfPitch,
                                                                    hexMesh.halfSide,
                                                                    order
                                                                )

    def getCoordsLines(self):
        return self.__coords, self.__lines

    def computeFirstOrderFunction(self, j):
        """
        Function to build first-order Wachspress elements (cf. Gout 1985).
        """
        line = lambda i: self.__lines['d%d' % hexCycle(i)]
        func = lambda x, y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y)
        a0, a1 = self.__coords['a%d' % j][0], self.__coords['a%d' % j][1]
        qai = computeCircleForRegularPolygon(self.__pitch, a0, a1)
        ci = qai / func(a0, a1)
        wi = lambda x, y: ci * func(x, y) / computeCircleForRegularPolygon(self.__pitch, x, y)

        return (wi,)

    def computeSecondOrderFunction(self, j):
        """
        Function to build second-order Wachspress elements (cf. Gout 1985).
        """
        line = lambda i: self.__lines['d%d' % hexCycle(i)]
        linep = self.__lines['dp%d' % j]
        func1 = lambda x, y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) * computeLi(linep, x, y)
        a0, a1 = self.__coords['a%d' % j][0], self.__coords['a%d' % j][1]
        qai = computeCircleForRegularPolygon(self.__pitch, a0, a1)
        ci = qai / func1(a0, a1)
        wi = lambda x, y: ci * func1(x, y) / computeCircleForRegularPolygon(self.__pitch, x, y)

        func2 = lambda x, y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) * computeLi(line(j), x, y)
        ip1 = lambda i: hexCycle(i)
        aiip1 = self.__coords['a%d%d' % (j, ip1(j+1))]
        qaiip1 = computeCircleForRegularPolygon(self.__pitch, aiip1[0], aiip1[1])
        ciip1 = qaiip1 / func2(aiip1[0], aiip1[1])
        wiip1 = lambda x, y: ciip1 * func2(x, y) / computeCircleForRegularPolygon(self.__pitch, x, y)

        return (wi, wiip1)

    def computeThirdOrderFunction(self, j):
        """
        Function to build third-order Wachspress elements (cf. Gout 1985).
        """
        line = lambda i: self.__lines['d%d' % hexCycle(i)]
        linep = lambda i, k: self.__lines['dp%d%d' % (hexCycle(i), hexCycle(k))]
        nextJ = hexCycle(j+1)
        func1 = lambda x, y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) * computeCircleForRegularPolygon(self.__pitch, x, y, 7.0/27.0)
        a0, a1 = self.__coords['a%d' % j][0], self.__coords['a%d' % j][1]
        qai = computeCircleForRegularPolygon(self.__pitch, a0, a1)
        ci = qai / func1(a0, a1)
        wi = lambda x, y: ci * func1(x, y) / computeCircleForRegularPolygon(self.__pitch, x, y)

        func2 = lambda x, y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) * computeLi(line(j), x, y) * computeLi(linep(j+1, j+3), x, y)
        aiiip1 = self.__coords['a%d%d%d' % (j, j, nextJ)]
        a0, a1 = aiiip1[0], aiiip1[1]
        qaiiip1 = computeCircleForRegularPolygon(self.__pitch, a0, a1)
        ciiip1 = qaiiip1 / func2(a0, a1)
        wiiip1 = lambda x, y: ciiip1 * func2(x, y) / computeCircleForRegularPolygon(self.__pitch, x, y)
        
        func3 = lambda x, y: computeLi(line(j+2), x, y) * computeLi(line(j+3), x, y) * computeLi(line(j+4), x, y) * computeLi(line(j+5), x, y) * computeLi(line(j), x, y) * computeLi(linep(j, j+4), x, y)
        aip1ii = self.__coords['a%d%d%d' % (nextJ,nextJ,j)]
        a0, a1 = aip1ii[0], aip1ii[1]
        qaip1ip1i = computeCircleForRegularPolygon(self.__pitch, a0, a1)
        cip1ip1i = qaip1ip1i / func3(a0, a1)
        wip1ip1i = lambda x, y: cip1ip1i * func3(x, y) / computeCircleForRegularPolygon(self.__pitch, x, y)

        return (wi, wiiip1, wip1ip1i)
