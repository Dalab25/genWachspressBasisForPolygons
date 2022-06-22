import numpy as np
import matplotlib.pyplot as plt

class IrregularPentagon(object):

    def __init__(self):
        self.__vertices = self.computeVertices()
        self.__coeffsAdjoint = self.computeAdjointIrregularPentagon()
        self.__dofsDict = None

    @property
    def vertices(self):
        return self.__vertices

    @property
    def coeffsAdjoint(self):
        return self.__coeffsAdjoint
    
    def computeVertices(self):    
        a_0 = [-2.5,-0.5]
        a_1 = [-1.5,2]
        a_2 = [2.5, 2]
        a_3 = [3.5,0.25]
        a_4 = [0,-2.5]
        vertices = [a_0, a_1, a_2, a_3, a_4]
        
        return vertices
    
    def getExternalIntersectionPoints(self):
        intersectionPoints = []
        intersectionPointL1L3 = [0.625/4.25, 2.5*(0.625/4.25) + 5.75]
        intersectionPointL2L4 = [15.75/2.75, 2]
        intersectionPointL3L5 = [8.875/0.95, -0.8*(8.875/0.95) - 2.5]
        intersectionPointL4L1 = [-28.875/6, -2.5*(28.875/6) + 5.75]
        intersectionPointL5L2 = [-5.625,2]
        
        intersectionPoints.append(intersectionPointL1L3)
        intersectionPoints.append(intersectionPointL2L4)
        intersectionPoints.append(intersectionPointL3L5)
        intersectionPoints.append(intersectionPointL4L1)
        intersectionPoints.append(intersectionPointL5L2)
        
        return intersectionPoints
    
    def computeAdjointIrregularPentagon(self):
        """
        Compute the adjoint, for that we consider a conic for the adjoint such as:
        x**2 + a*y**2 + b*x*y + c*x + d*y + e
        Consequently, for a given list of points (x_i,y_i), conic(x_i,y_i) = 0 
        <=> a*y_i**2 + b*x_i*y_i + c*x_i + d*y_i + e = -x_i**2
        
        Hence, we can create a matrix A such as:
        A = [y_i**2, x_i*y_i, x_i, y_i, 1]
        
        Consequently, a linear system can be written with the following vectors
        X = [a,b,c,d,e]
        b = [-xi**2]
        
        AX = b <=> X = A^-1b
        """
        intersectionPoints = self.getExternalIntersectionPoints()
        A = np.zeros((5,5))
        b = np.zeros(5)
        
        for i in range(len(intersectionPoints)):
            A[i][0] = intersectionPoints[i][1]**2
            A[i][1] = intersectionPoints[i][0]*intersectionPoints[i][1]
            A[i][2] = intersectionPoints[i][0]
            A[i][3] = intersectionPoints[i][1]
            A[i][4] = 1
            b[i] = -intersectionPoints[i][0]**2 #-x_i**2
            
        coeffs = np.linalg.solve(A, b)
        return coeffs
    
      
    def plotMeshFromVertices(self):
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        vertices = [v for v in self.__vertices]
        vertices.append(vertices[0])
        plt.plot(*zip(*vertices))
        
        intersectionPoints = self.getExternalIntersectionPoints()
        x = []
        y = []
        
        for i in range(len(intersectionPoints)):
            x.append(intersectionPoints[i][0])
            y.append(intersectionPoints[i][1])
        plt.scatter(x,y)
            
        #Tracer la premiere conique
        x = np.linspace(-10, 10, 800)
        y = np.linspace(-20, 15, 800)
        x, y = np.meshgrid(x, y)
        
        coeffs = self.coeffsAdjoint
        a,b,c,d,e = coeffs[0], coeffs[1], coeffs[2], coeffs[3],coeffs[4]
        print(a,b,c,d,e)
        plt.contour(x, y,(x**2 + a*y**2 + b*x*y +  c*x + d*y + e), [0], colors='r')
            
        xs = []
        ys = []
        if self.__dofsDict:
            for i in range(self.__nSides):
                coordsOnSide = self.__dofsDict["side{}".format(i)]
                xtmp, ytmp = zip(*coordsOnSide)
                xs.extend(xtmp)
                ys.extend(ytmp)
            plt.scatter(xs,ys)
                
        plt.tight_layout()
        plt.show() 
