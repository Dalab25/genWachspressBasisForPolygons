import numpy as np

import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib import cm

class PolygonForWachspress(object):

    def __init__(self, nSides, radius):
        self.__nSides = nSides
        self.__radius = radius # 2hs in solver
        self.__deltaAngle = 2*np.pi/self.__nSides

        self.__vertices = self.computePolygonVertices()
        self.__halfSide = self.__radius * np.sin(self.__deltaAngle/2)
        self.__sidelength = 2 * self.__halfSide
        #self.__halfSide = 0.5 * self.__sidelength 
        #self.__sidelength = abs(self.__vertices[0][0])*2
        #self.__halfSide = 0.5 * self.__sidelength / (3**0.5)
        self.__dofsDict = None

    @property
    def nSides(self):
        return self.__nSides
    
    @property
    def vertices(self):
        return self.__vertices

    @property
    def pitch(self):
        return 2*self.halfPitch
    
    @property
    def side(self):
        return self.__sidelength

    @property
    def halfPitch(self):
        return self.__radius * np.cos(self.__deltaAngle/2) if self.__nSides % 2 == 0 else None

    @property
    def adjointRadius(self):
        """
        The radius of the adjoint polygon, always a circle for regular
        polygons.  It is equivalent to twice the distance of the
        centre of the polygon to a given side of the polygon (for
        polygons with even number of sides, it is equivalent to the
        pitch."
        """
        return self.__radius * np.cos(self.__deltaAngle/2) * 2

    @property
    def halfSide(self):
        return self.__halfSide

    def setDofs(self, dofsDict):
        self.__dofsDict = dofsDict 
    
    def computePolygonVertices(self):
        xc = 0.
        yc = 0.
        deltaAngle = 2 * np.pi / self.__nSides
        # start angle
        angle = (np.pi - deltaAngle)/2 + np.pi

        vertices = []

        for i in range(self.__nSides):
            x = xc + self.__radius * np.cos(angle)
            y = yc + self.__radius * np.sin(angle)
            angle += deltaAngle
            vertices.append((x,y))

        return vertices

    def plotPolygon(self):
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        patches = []
        pt = self.__vertices
        hexagon = Polygon(pt, True)
        patches.append(hexagon)
        ax.annotate(len(patches)-1, ((pt[0][0]+pt[3][0])/2, (pt[0][1]+pt[3][1])/2), color='black', weight='bold', fontsize=6, ha='center', va='center')

        p = PatchCollection(patches, cmap=cm.jet)
        colors = 100*np.random.rand(len(patches))
        p.set_array(np.array(colors))
        ax.add_collection(p)
        ax.axis('scaled')  # ax.autoscale(enable=True)
        plt.colorbar(p)

        plt.tight_layout()
        plt.show()

    def plotMeshFromVertices(self):
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        # mpl.rcParams['lines.color'] = 'k'
        # mpl.rcParams['axes.prop_cycle'] = mpl.cycler('color', ['k'])
        # nx = ny = 1.1
        # xWindow = -nx*side, nx*side
        # yWindow = -ny*hp, ny*hp
        # N = 400
        # xPoints = np.linspace(*xWindow, N)
        # yPoints = np.linspace(*yWindow, N)
        # x, y = np.meshgrid(xPoints, yPoints)
        # axes()   
        # # ax.axis('scaled')
        vertices = [v for v in self.__vertices]
        vertices.append(vertices[0])
        plt.plot(*zip(*vertices))
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


def main():
    print("Running")
    n = 10
    radius = 1.
    poly = PolygonForWachspress(n, radius)
    order = 4

    vertices = poly.vertices
    print(vertices)
    poly.plotMeshFromVertices()
    #poly.plotPolygon()
    """
    vDict = {}
    for i in range(n):
        verticesOnSide = [vertices[i]]
        a = vertices[i]
        b = vertices[nCycle(i+1, n-1)]
        for j in range(1, order):
            coords = computeBarycentre2(a, b, order, j)
            verticesOnSide.append(coords)
        vDict["side{}".format(i)] = verticesOnSide
    for key, value in vDict.items():
        print(key, value)
    poly.setDofs(vDict)
    #poly.plotMeshFromVertices()
    """

if __name__ == "__main__":
    main()
