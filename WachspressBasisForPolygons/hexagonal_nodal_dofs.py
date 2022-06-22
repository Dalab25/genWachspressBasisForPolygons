# -*- coding: utf-8 -*-
#
# File: hexagonal_nodal_dofs.py
#
import numpy as np
from geometrical_elements import computeLineCoefficients, computeBarycentre, computeBarycentre2

def buildCoordsAndLinesForWSF(hp, hs, order):
    """
    Evaluates the geometrical constructs (nodes and line segments)
    to build the Wachspress basis functions up to order 5.
    
    It should respect the following order for a given order k:
    List[a_i], List[a_i_1], List[a_i_2], List[a_i..._k-1]
    Example at the order 3:
    a1,a2,a3..a6; a_1_1,a_2_1,..a_6_1; a_1_2,a_2_2, .. a_6_2.
    
    Witout this order, the calculations will be impossible.
    """
    a1 = (-hs, -hp)
    a2 = (hs, -hp)
    a3 = (2*hs, 0.)
    a4 = (hs, hp)
    a5 = (-hs, hp)
    a6 = (-2*hs, 0.)

    # order 1
    coordsDict = {'a1': a1,
                  'a2': a2,
                  'a3': a3,
                  'a4': a4,
                  'a5': a5,
                  'a6': a6,
                 }

    d1 = computeLineCoefficients(a1, a6)
    d2 = computeLineCoefficients(a2, a1)
    d3 = computeLineCoefficients(a3, a2)
    d4 = computeLineCoefficients(a4, a3)
    d5 = computeLineCoefficients(a5, a4)
    d6 = computeLineCoefficients(a6, a5)

    linesDict = {'d1': d1,
                 'd2': d2,
                 'd3': d3,
                 'd4': d4,
                 'd5': d5,
                 'd6': d6,
                 }

    # order 2
    if order == 2:
        coordsDict['a12'] = computeBarycentre(a1, a2, 1., 1.)
        coordsDict['a23'] = computeBarycentre(a2, a3, 1., 1.)
        coordsDict['a34'] = computeBarycentre(a3, a4, 1., 1.)
        coordsDict['a45'] = computeBarycentre(a4, a5, 1., 1.)
        coordsDict['a56'] = computeBarycentre(a5, a6, 1., 1.)
        coordsDict['a61'] = computeBarycentre(a6, a1, 1., 1.)

        linesDict['dp1'] = computeLineCoefficients(coordsDict['a61'], coordsDict['a12'])
        linesDict['dp2'] = computeLineCoefficients(coordsDict['a12'], coordsDict['a23'])
        linesDict['dp3'] = computeLineCoefficients(coordsDict['a23'], coordsDict['a34'])
        linesDict['dp4'] = computeLineCoefficients(coordsDict['a34'], coordsDict['a45'])
        linesDict['dp5'] = computeLineCoefficients(coordsDict['a45'], coordsDict['a56'])
        linesDict['dp6'] = computeLineCoefficients(coordsDict['a56'], coordsDict['a61'])

    # order 3
    if order == 3:
        """
        coordsDict['a112'] = computeBarycentre(a1, a2, 2., 1.)
        coordsDict['a221'] = computeBarycentre(a2, a1, 2., 1.)
        coordsDict['a223'] = computeBarycentre(a2, a3, 2., 1.)
        coordsDict['a332'] = computeBarycentre(a3, a2, 2., 1.)
        coordsDict['a334'] = computeBarycentre(a3, a4, 2., 1.)
        coordsDict['a443'] = computeBarycentre(a4, a3, 2., 1.)
        coordsDict['a445'] = computeBarycentre(a4, a5, 2., 1.)
        coordsDict['a554'] = computeBarycentre(a5, a4, 2., 1.)
        coordsDict['a556'] = computeBarycentre(a5, a6, 2., 1.)
        coordsDict['a665'] = computeBarycentre(a6, a5, 2., 1.)
        coordsDict['a661'] = computeBarycentre(a6, a1, 2., 1.)
        coordsDict['a116'] = computeBarycentre(a1, a6, 2., 1.)
        """
        #"""
        coordsDict['a112'] = computeBarycentre(a1, a2, 2., 1.)
        coordsDict['a223'] = computeBarycentre(a2, a3, 2., 1.)
        coordsDict['a334'] = computeBarycentre(a3, a4, 2., 1.)
        coordsDict['a445'] = computeBarycentre(a4, a5, 2., 1.)
        coordsDict['a556'] = computeBarycentre(a5, a6, 2., 1.)
        coordsDict['a661'] = computeBarycentre(a6, a1, 2., 1.)
        coordsDict['a221'] = computeBarycentre(a2, a1, 2., 1.)
        coordsDict['a332'] = computeBarycentre(a3, a2, 2., 1.)
        coordsDict['a443'] = computeBarycentre(a4, a3, 2., 1.)
        coordsDict['a554'] = computeBarycentre(a5, a4, 2., 1.)
        coordsDict['a665'] = computeBarycentre(a6, a5, 2., 1.)
        coordsDict['a116'] = computeBarycentre(a1, a6, 2., 1.)
        #"""
        linesDict['dp15'] = computeLineCoefficients(coordsDict['a112'], coordsDict['a554'])
        linesDict['dp24'] = computeLineCoefficients(coordsDict['a221'], coordsDict['a445'])
        linesDict['dp26'] = computeLineCoefficients(coordsDict['a223'], coordsDict['a665'])
        linesDict['dp35'] = computeLineCoefficients(coordsDict['a332'], coordsDict['a556'])
        linesDict['dp31'] = computeLineCoefficients(coordsDict['a334'], coordsDict['a116'])
        linesDict['dp46'] = computeLineCoefficients(coordsDict['a443'], coordsDict['a661'])
        linesDict['dp51'] = linesDict['dp15']
        linesDict['dp42'] = linesDict['dp24']
        linesDict['dp62'] = linesDict['dp26']
        linesDict['dp53'] = linesDict['dp35']
        linesDict['dp13'] = linesDict['dp31']
        linesDict['dp64'] = linesDict['dp46']

    # order 4
    if order == 4:
        """
        coordsDict['a112'] = computeBarycentre(a1, a2, 3., 1.)
        coordsDict['a221'] = computeBarycentre(a2, a1, 3., 1.)
        coordsDict['a223'] = computeBarycentre(a2, a3, 3., 1.)
        coordsDict['a332'] = computeBarycentre(a3, a2, 3., 1.)
        coordsDict['a334'] = computeBarycentre(a3, a4, 3., 1.)
        coordsDict['a443'] = computeBarycentre(a4, a3, 3., 1.)
        coordsDict['a445'] = computeBarycentre(a4, a5, 3., 1.)
        coordsDict['a554'] = computeBarycentre(a5, a4, 3., 1.)
        coordsDict['a556'] = computeBarycentre(a5, a6, 3., 1.)
        coordsDict['a665'] = computeBarycentre(a6, a5, 3., 1.)
        coordsDict['a661'] = computeBarycentre(a6, a1, 3., 1.)
        coordsDict['a116'] = computeBarycentre(a1, a6, 3., 1.)

        coordsDict['a12'] = computeBarycentre(a1, a2, 1., 1.)
        coordsDict['a23'] = computeBarycentre(a2, a3, 1., 1.)
        coordsDict['a34'] = computeBarycentre(a3, a4, 1., 1.)
        coordsDict['a45'] = computeBarycentre(a4, a5, 1., 1.)
        coordsDict['a56'] = computeBarycentre(a5, a6, 1., 1.)
        coordsDict['a61'] = computeBarycentre(a6, a1, 1., 1.)
        """
        
        coordsDict['a112'] = computeBarycentre(a1, a2, 3., 1.)
        coordsDict['a223'] = computeBarycentre(a2, a3, 3., 1.)
        coordsDict['a334'] = computeBarycentre(a3, a4, 3., 1.)
        coordsDict['a445'] = computeBarycentre(a4, a5, 3., 1.)
        coordsDict['a556'] = computeBarycentre(a5, a6, 3., 1.)
        coordsDict['a661'] = computeBarycentre(a6, a1, 3., 1.)
        coordsDict['a12'] = computeBarycentre(a1, a2, 1., 1.)
        coordsDict['a23'] = computeBarycentre(a2, a3, 1., 1.)
        coordsDict['a34'] = computeBarycentre(a3, a4, 1., 1.)
        coordsDict['a45'] = computeBarycentre(a4, a5, 1., 1.)
        coordsDict['a56'] = computeBarycentre(a5, a6, 1., 1.)
        coordsDict['a61'] = computeBarycentre(a6, a1, 1., 1.)
        coordsDict['a221'] = computeBarycentre(a2, a1, 3., 1.)
        coordsDict['a332'] = computeBarycentre(a3, a2, 3., 1.)
        coordsDict['a443'] = computeBarycentre(a4, a3, 3., 1.)
        coordsDict['a554'] = computeBarycentre(a5, a4, 3., 1.)
        coordsDict['a665'] = computeBarycentre(a6, a5, 3., 1.)
        coordsDict['a116'] = computeBarycentre(a1, a6, 3., 1.)
        
    # order 5
    if order > 4:
        """
        coordsDict['a11'] = computeBarycentre2(a1, a2, 5, 1)
        coordsDict['a12'] = computeBarycentre2(a1, a2, 5, 2)
        coordsDict['a13'] = computeBarycentre2(a1, a2, 5, 3)
        coordsDict['a14'] = computeBarycentre2(a1, a2, 5, 4)

        coordsDict['a21'] = computeBarycentre2(a2, a3, 5, 1)
        coordsDict['a22'] = computeBarycentre2(a2, a3, 5, 2)
        coordsDict['a23'] = computeBarycentre2(a2, a3, 5, 3)
        coordsDict['a24'] = computeBarycentre2(a2, a3, 5, 4)

        coordsDict['a31'] = computeBarycentre2(a3, a4, 5, 1)
        coordsDict['a32'] = computeBarycentre2(a3, a4, 5, 2)
        coordsDict['a33'] = computeBarycentre2(a3, a4, 5, 3)
        coordsDict['a34'] = computeBarycentre2(a3, a4, 5, 4)

        coordsDict['a41'] = computeBarycentre2(a4, a5, 5, 1)
        coordsDict['a42'] = computeBarycentre2(a4, a5, 5, 2)
        coordsDict['a43'] = computeBarycentre2(a4, a5, 5, 3)
        coordsDict['a44'] = computeBarycentre2(a4, a5, 5, 4)

        coordsDict['a51'] = computeBarycentre2(a5, a6, 5, 1)
        coordsDict['a52'] = computeBarycentre2(a5, a6, 5, 2)
        coordsDict['a53'] = computeBarycentre2(a5, a6, 5, 3)
        coordsDict['a54'] = computeBarycentre2(a5, a6, 5, 4)

        coordsDict['a61'] = computeBarycentre2(a6, a1, 5, 1)
        coordsDict['a62'] = computeBarycentre2(a6, a1, 5, 2)
        coordsDict['a63'] = computeBarycentre2(a6, a1, 5, 3)
        coordsDict['a64'] = computeBarycentre2(a6, a1, 5, 4)
        """
        coordsDict['a11'] = computeBarycentre2(a1, a2, 5, 1)
        coordsDict['a21'] = computeBarycentre2(a2, a3, 5, 1)
        coordsDict['a31'] = computeBarycentre2(a3, a4, 5, 1)
        coordsDict['a41'] = computeBarycentre2(a4, a5, 5, 1)
        coordsDict['a51'] = computeBarycentre2(a5, a6, 5, 1)
        coordsDict['a61'] = computeBarycentre2(a6, a1, 5, 1)
        

        coordsDict['a12'] = computeBarycentre2(a1, a2, 5, 2)
        coordsDict['a22'] = computeBarycentre2(a2, a3, 5, 2)
        coordsDict['a32'] = computeBarycentre2(a3, a4, 5, 2)
        coordsDict['a42'] = computeBarycentre2(a4, a5, 5, 2)
        coordsDict['a52'] = computeBarycentre2(a5, a6, 5, 2)
        coordsDict['a62'] = computeBarycentre2(a6, a1, 5, 2)
        
        coordsDict['a13'] = computeBarycentre2(a1, a2, 5, 3)
        coordsDict['a23'] = computeBarycentre2(a2, a3, 5, 3)
        coordsDict['a33'] = computeBarycentre2(a3, a4, 5, 3)
        coordsDict['a43'] = computeBarycentre2(a4, a5, 5, 3)
        coordsDict['a53'] = computeBarycentre2(a5, a6, 5, 3)
        coordsDict['a63'] = computeBarycentre2(a6, a1, 5, 3)

        coordsDict['a14'] = computeBarycentre2(a1, a2, 5, 4)
        coordsDict['a24'] = computeBarycentre2(a2, a3, 5, 4)
        coordsDict['a34'] = computeBarycentre2(a3, a4, 5, 4)
        coordsDict['a44'] = computeBarycentre2(a4, a5, 5, 4)
        coordsDict['a54'] = computeBarycentre2(a5, a6, 5, 4)
        coordsDict['a64'] = computeBarycentre2(a6, a1, 5, 4)
    return coordsDict, linesDict
