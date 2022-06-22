# -*- coding: utf-8 -*-
import unittest
import numpy as np
from WSF_to_generate_irregular_pentagon import WSFTogenerateIrregularPentagon
from irregular_pentagon import IrregularPentagon
from geometrical_elements import pentaCycle

class WachspressShapeFunctionsTests(unittest.TestCase):
    def test_computeSecondOrderIrregularPentagon(self):
        print("Running")
        order = 2
        pentagon = IrregularPentagon()
        wsf1 = WSFTogenerateIrregularPentagon(order, pentagon)
        func = wsf1.computeSecondOrderFunction

        #vertices = wsf1.pentagon.vertices
        #verticesOnSide = wsf1.verticesOnSide
        coords = wsf1.coords
        eps = 10**-14
        print("Verification of the equality wi(aj) = delta_ij, wi(aj_1) = 0")
        print("And also, wiip1(aj) = 0, wiip1(aj_1) = delta_ij")
        ip1 = lambda i: pentaCycle(i)
        for i in range(5):
            (wi,wiip1) = func(i)
            for j in range(5):
                ajjx, ajjy = coords['a%d%d' % (j, ip1(j+1))][0], coords['a%d%d' % (j, ip1(j+1))][1]
                ajx,ajy = coords['a%d' % j][0], coords['a%d' % j][1]
                if(i==j):
                    assert(np.abs(wi(ajx, ajy) - 1) < eps)
                    assert(np.abs(wiip1(ajjx,ajjy)) - 1 < eps)
                else:
                    assert(wi(ajx,ajy) == 0)
                    assert(wi(ajjx,ajjy) == 0)
                    assert(wiip1(ajx,ajy) == 0)
                    assert(wiip1(ajjx,ajjy) == 0)

if __name__ == '__main__':
    unittest.main()
