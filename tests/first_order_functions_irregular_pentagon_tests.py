# -*- coding: utf-8 -*-
import unittest
from WSF_to_generate_irregular_pentagon import WSFTogenerateIrregularPentagon
from irregular_pentagon import IrregularPentagon
import numpy as np


class WachspressShapeFunctionsTests(unittest.TestCase):
    def test_computeFirstOrderIrregularPentagon(self):
        print("Running")
        order = 1

        #poly = PolygonForWachspress(n, -0.8660254037844388)
        pentagon = IrregularPentagon()
        wsf1 = WSFTogenerateIrregularPentagon(order, pentagon)
        func = wsf1.computeFirstOrderFunction

        vertices = wsf1.pentagon.vertices
        print(vertices)
        coords = wsf1.coords
        #print(vertices)


        eps = 10**-14
        print("Test 1) verification of the equality wi(ai) = 1, wi(aj) = 0 for a hexagon")
        for i in range(5):
            wi = func(i)
            #wi2 = func(i + 1, poly.pitch, vertices, lines)
            for j in range(5):
                ajx,ajy = vertices[j][0],vertices[j][1]
                print("wi(ajx, ajy): ", wi(ajx,ajy))
                if(i == j):
                    assert(np.abs(wi(ajx,ajy) -1)  < eps)
                else:
                    assert(wi(ajx,ajy) == 0)
                #print("Hoop:", wi2(ajx,ajy))
            print("")
            #assert(wi(aix, aiy) - 1 < eps)

if __name__ == '__main__':
    unittest.main()
