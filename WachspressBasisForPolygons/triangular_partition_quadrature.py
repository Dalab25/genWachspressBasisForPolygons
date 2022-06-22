# -*- coding: utf-8 -*-
#
# File: triangular_partition_quadrature.py
#
"""
Taken from Rathod 2007: "On the application of two symmetric
Gauss-Legendre quadrature rules for composite numerical integration
over a triangular surface", Applied Mathematics and Computation 190
(2007)
"""
import numpy as np


def transformPoint(r, s, tri):
    """
    R:
    (0,1) |
          |\
          | \
          |  \
    (0,0) |---(1,0)-
    Tranform from reference triangle R to physical triangle T:
    |x| = F(r,s) = | ax - cx   bx - cx | |r| + |cx|
    |y|            | ay - cy   by - cy | |s| + |cy|
    where (ax,ay), (bx,by) and (cx,cy) are the nodes of T
    """
    x = (tri[0][0] - tri[2][0]) * r + (tri[1][0] - tri[2][0]) * s + tri[2][0]
    y = (tri[0][1] - tri[2][1]) * r + (tri[1][1] - tri[2][1]) * s + tri[2][1]
    return (x, y)


class TriangularPartitionQuadrature(object):

    @staticmethod
    def getNodesAndWeights(order, pitch):
        mi, wi = np.polynomial.legendre.leggauss(order)
        mj, wj = np.polynomial.legendre.leggauss(order)

        hp = 0.5 * pitch
        hs = 0.5 * pitch / (3**0.5)
        triangles = [[(0, 0), (-hs, -hp), (hs, -hp)],
                     [(0, 0), (hs, -hp), (2*hs, 0)],
                     [(0, 0), (2*hs, 0), (hs, hp)],
                     [(0, 0), (hs, hp), (-hs, hp)],
                     [(0, 0), (-hs, hp), (-2*hs, 0)],
                     [(0, 0), (-2*hs, 0), (-hs, -hp)]
                     ]

        hexArea = 3**0.5 * 0.5 * pitch**2  # area of regular hexagon
        area = 0.5 * hp * (hs*2)  # area of equilateral triangle in regular hexagon
        facw = area / 0.5

        quadrature = np.empty(6 * order * order * 3, dtype=object).reshape(6 * order * order, 3)
        weights = 0.
        iq = 0
        for tri in triangles:
            for i in range(order):
                for j in range(order):
                    wk = (1+mi[i]) * wi[i] * wj[j] / 8
                    wk *= facw
                    weights += wk
                    xk = (1+mi[i]) * (1+mj[j]) / 4
                    yk = (1+mi[i]) * (1-mj[j]) / 4
                    xt, yt = transformPoint(xk, yk, tri)
                    quadrature[iq] = (xt, yt, wk)
                    iq += 1

        assert(abs(hexArea/weights - 1.) < 1.e-14)

        return quadrature.reshape(len(quadrature), 3)
