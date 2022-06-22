# -*- coding: utf-8 -*-
#
# File: pentagonal_nodal_dofs.py
#
from geometrical_elements import computeLineCoefficients, computeBarycentre

def buildCoordsAndLinesForWSF(irregularPentagon, order):
    """
    Evaluates the geometrical constructs (nodes and line segments)
    to build the Wachspress basis functions up to order 4
    (cf. Gout1985)
    """
    vertices = irregularPentagon.vertices

    a0 = vertices[0]
    a1 = vertices[1]
    a2 = vertices[2]
    a3 = vertices[3]
    a4 = vertices[4]

    # order 1
    coordsDict = {'a0': a0,
                    'a1': a1,
                    'a2': a2,
                    'a3': a3,
                    'a4': a4
    }

    d0 = computeLineCoefficients(a0, a4, verbose=0)
    d1 = computeLineCoefficients(a1, a0, verbose=0)
    d2 = computeLineCoefficients(a2, a1, verbose=0)
    d3 = computeLineCoefficients(a3, a2, verbose=0)
    d4 = computeLineCoefficients(a4, a3, verbose=0)

    linesDict = {'d0': d0,
                    'd1': d1,
                    'd2': d2,
                    'd3': d3,
                    'd4': d4
    }

    # order 2
    if order > 1:
        coordsDict['a01'] = computeBarycentre(a0, a1, 1., 1.)
        coordsDict['a12'] = computeBarycentre(a1, a2, 1., 1.)
        coordsDict['a23'] = computeBarycentre(a2, a3, 1., 1.)
        coordsDict['a34'] = computeBarycentre(a3, a4, 1., 1.)
        coordsDict['a40'] = computeBarycentre(a4, a0, 1., 1.)

        linesDict['dp0'] = computeLineCoefficients(coordsDict['a40'], coordsDict['a01'])
        linesDict['dp1'] = computeLineCoefficients(coordsDict['a01'], coordsDict['a12'])
        linesDict['dp2'] = computeLineCoefficients(coordsDict['a12'], coordsDict['a23'])
        linesDict['dp3'] = computeLineCoefficients(coordsDict['a23'], coordsDict['a34'])
        linesDict['dp4'] = computeLineCoefficients(coordsDict['a34'], coordsDict['a40'])

    # order 3
    if order > 2:
        #'''
        coordsDict['a001'] = computeBarycentre(a0, a1, 2., 1.)
        coordsDict['a110'] = computeBarycentre(a1, a0, 2., 1.)
        coordsDict['a112'] = computeBarycentre(a1, a2, 2., 1.)
        coordsDict['a221'] = computeBarycentre(a2, a1, 2., 1.)
        coordsDict['a223'] = computeBarycentre(a2, a3, 2., 1.)
        coordsDict['a332'] = computeBarycentre(a3, a2, 2., 1.)
        coordsDict['a334'] = computeBarycentre(a3, a4, 2., 1.)
        coordsDict['a443'] = computeBarycentre(a4, a3, 2., 1.)
        coordsDict['a440'] = computeBarycentre(a4, a0, 2., 1.)
        coordsDict['a004'] = computeBarycentre(a0, a4, 2., 1.)

        linesDict['dp04'] = computeLineCoefficients(coordsDict['a001'], coordsDict['a443'])
        linesDict['dp13'] = computeLineCoefficients(coordsDict['a110'], coordsDict['a334'])
        linesDict['dp14'] = computeLineCoefficients(coordsDict['a112'], coordsDict['a004'])
        linesDict['dp24'] = computeLineCoefficients(coordsDict['a221'], coordsDict['a440'])
        linesDict['dp20'] = computeLineCoefficients(coordsDict['a223'], coordsDict['a004'])
        linesDict['dp34'] = computeLineCoefficients(coordsDict['a332'], coordsDict['a440'])
        linesDict['dp40'] = linesDict['dp04']
        linesDict['dp31'] = linesDict['dp13']
        linesDict['dp41'] = linesDict['dp14']
        linesDict['dp42'] = linesDict['dp24']
        linesDict['dp02'] = linesDict['dp20']
        linesDict['dp43'] = linesDict['dp34']
        #'''
    #order 4
    if order > 3:
        coordsDict['a001'] = computeBarycentre(a0, a1, 3., 1.)
        coordsDict['a110'] = computeBarycentre(a1, a0, 3., 1.)
        coordsDict['a112'] = computeBarycentre(a1, a2, 3., 1.)
        coordsDict['a221'] = computeBarycentre(a2, a1, 3., 1.)
        coordsDict['a223'] = computeBarycentre(a2, a3, 3., 1.)
        coordsDict['a332'] = computeBarycentre(a3, a2, 3., 1.)
        coordsDict['a334'] = computeBarycentre(a3, a4, 3., 1.)
        coordsDict['a443'] = computeBarycentre(a4, a3, 3., 1.)
        coordsDict['a440'] = computeBarycentre(a4, a0, 3., 1.)
        coordsDict['a004'] = computeBarycentre(a0, a4, 3., 1.)

        coordsDict['a01'] = computeBarycentre(a0, a1, 1., 1.)
        coordsDict['a12'] = computeBarycentre(a1, a2, 1., 1.)
        coordsDict['a23'] = computeBarycentre(a2, a3, 1., 1.)
        coordsDict['a34'] = computeBarycentre(a3, a4, 1., 1.)
        coordsDict['a40'] = computeBarycentre(a4, a0, 1., 1.)
        
    return coordsDict, linesDict