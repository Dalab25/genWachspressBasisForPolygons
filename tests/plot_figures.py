# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from irregular_pentagon import IrregularPentagon

def plot_order4(boolConstraints):
    
    if(boolConstraints):
        solution = np.loadtxt("../results/polygon_5sides/order4/symbolic_coefficients/solution_order_4_constraints1.txt")
    else:
        solution = np.loadtxt("../results/polygon_5sides/order4/symbolic_coefficients/solution_order_4_no_constraints1.txt")      
    mesh = IrregularPentagon()
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    vertices = [v for v in mesh.computeVertices()]
    vertices.append(vertices[0])
    plt.plot(*zip(*vertices))
    #Tracer la premiere conique
    x = np.linspace(-10, 10, 800)
    y = np.linspace(-20, 15, 800)
    x, y = np.meshgrid(x, y)
    
    #a,b,c,d,e = coeffs[0], coeffs[1], coeffs[2], coeffs[3],coeffs[4]
    a,b,c,d,e,f,g,h,i,j = solution[0:10]
    plt.contour(x, y,(a*x**3 + b*y**3 + c*(x**2)*y + d*x*(y**2) + e*x**2 + f*y**2 + g*x*y + h*x + i*y + j), [0], colors='r')
    
    plt.xlim(-2.55,0.05)
    plt.ylim(-3,3)
    #ax.set_yticklabels(y, fontsize=5)
    #ax.set_xticklabels(x, fontsize=5)
    #plt.tight_layout()
    plt.show() 
    
def plot_order3(boolConstraints):
    
    if(boolConstraints):
        solution = np.loadtxt("../results/polygon_5sides/order3/symbolic_coefficients/solution_order_3_constraints1.txt")
    else:
        solution = np.loadtxt("../results/polygon_5sides/order3/symbolic_coefficients/solution_order_3_no_constraints1.txt")      

    mesh = IrregularPentagon()
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    vertices = [v for v in mesh.computeVertices()]
    vertices.append(vertices[0])
    plt.plot(*zip(*vertices))
    #Tracer la premiere conique
    x = np.linspace(-10, 10, 800)
    y = np.linspace(-20, 15, 800)
    x, y = np.meshgrid(x, y)
    
    #a,b,c,d,e = coeffs[0], coeffs[1], coeffs[2], coeffs[3],coeffs[4]
    a,b,c,d,e,f = solution[0:6]
    plt.contour(x, y,(a*x**2 + b*y**2 + c*x*y +  d*x + e*y + f), [0], colors='r')

    plt.xlim(-2.55,0.05)
    plt.ylim(-3,3)
    #plt.tight_layout()
    plt.show() 

def main():
    bool_constraints = False
    plot_order3(bool_constraints)
    #plot_order4(bool_constraints)
    #plot_order3(bool_constraints)
    
if __name__ == "__main__":
    main()
