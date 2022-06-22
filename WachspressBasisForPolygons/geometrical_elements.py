# -*- coding: utf-8 -*-
#
# File: geometrical_elements.py
#

def computeLineCoefficients(p, q, verbose=0):
    """
    params: p and q are tuples (2D) of the form (x,y).
    line joining p and q is of form ax + by + c = 0
    """
    a = p[1]-q[1]
    b = -(p[0]-q[0])
    c = p[0]*q[1] - q[0]*p[1]
    if verbose:
        print("Line equation is ", a, "x + ", b, "y + ", c, " = 0")
    return a, b, c

def computeCircleForRegularPolygon(p, x, y, fac=1.):
    """
    Evaluates the circle centred on (0,0) of radius p and a given
    scaling factor fac.
    """
    return x**2 + y**2 - fac*p**2

def computeCircleForRegularHex(p, x, y, fac=1.):
    """
    Evaluates the circle centred on (0,0) of radius p and a given
    scaling factor fac.
    """
    return x**2 + y**2 - fac*p**2


def computeAdjointForIrregularPentagon(coeffs,x,y):
    """
    Evaluates the adjoint of the irregular convex pentagon
    """
    return x**2 + coeffs[0]*y**2 + coeffs[1]*x*y + coeffs[2]*x + coeffs[3]*y + coeffs[4]


def nCycle(i, n):
    return n-1 if  (i+1)%n == 0 else (i+1)%n - 1

def computeBarycentre2(p, q, order, j):
    """
    Computes weighted centre of p and q for a given order and for
    nodes in increasing j from p.
    """
    a_i_j_x = (1/order)*((order - j)*p[0] + j*q[0])
    a_i_j_y = (1/order)*((order - j)*p[1] + j*q[1])
    return (a_i_j_x, a_i_j_y)

def computeBarycentre(p, q, wp, wq):
    wpn = wp/(wp+wq)
    wqn = wq/(wp+wq)
    return ((wpn*p[0]+wqn*q[0]), (wpn*p[1]+wqn*q[1]))

def computeFunctionOrderK(coeffs,order,x,y):
    sumTerms = 0
    compteur = 0
    degrees = [(i, j) for i in range(order+1) for j in range(order+1) if i+j <= order]
    for beta, gamma in degrees:
        sumTerms += coeffs[compteur]*(x**beta)*(y**gamma)
        compteur += 1
    return sumTerms

 
def computeCircleForRegularHexCoeffs(coeffs, x, y, COEFFS3=True):
    """
    Evaluates the circle centred on (0,0) of radius p and a given
    scaling factor fac.
    """
    if COEFFS3:
        return coeffs[0]*x**2 + coeffs[1]*y**2 + coeffs[2]
    else:
        return coeffs[0]*x**2 + coeffs[1]*y**2 + coeffs[2]*x*y + coeffs[3]*x + coeffs[4]*y + coeffs[5]

def pentaCycle(i):
    """
    Returns the corresponding edge for an irregular pentagon.
    """
    return 4 if (i+1) % 5 == 0 else (i+1)%5 - 1

def computeLi(coeffs, x, y):
    """
    Evaluates the line equation ax + by + c
    """
    return coeffs[0]*x + coeffs[1]*y + coeffs[2]


def computeQuartic(coeffs, x, y):
    """
    Evaluates the quartic polynomial in 2D
    """
    return coeffs[0]*x**4 + coeffs[1]*y**4 + coeffs[2]*(x**3)*y + coeffs[3]*x*(y**3) + coeffs[4]*(x**2)*y**2 + coeffs[5]*x**3 + coeffs[6]*y**3 + coeffs[7]*(x**2)*y + coeffs[8]*x*(y**2) + coeffs[9]*x**2 + coeffs[10]*y**2 + coeffs[11]*x*y + coeffs[12]*x + coeffs[13]*y + coeffs[14]


def computeCubic(coeffs, x, y):
    """
    Evaluates the cubic equation in 2D
    """
    return coeffs[0]*x**3 + coeffs[1]*y**3 + coeffs[2]*(x**2)*y + coeffs[3]*x*(y**2) + coeffs[4]*x**2 + coeffs[5]*y**2 + coeffs[6]*x*y + coeffs[7]*x + coeffs[8]*y + coeffs[9]


def computeConic(coeffs, x, y):
    """
    Evaluates the conic equation
    """
    return coeffs[0]*x**2 + coeffs[1]*y**2 + coeffs[2]*x*y + coeffs[3]*x + coeffs[4]*y + coeffs[5]

def hexCycle(i):
    """
    Returns the corresponding edge for a regular hexagon.
    """
    return 6 if i % 6 == 0 else (i % 6)