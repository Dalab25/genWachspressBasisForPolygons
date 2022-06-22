# -*- coding: utf-8 -*-
#
# File: honeycomb.py
#
import warnings
import numpy as np

class HoneyCombMesh(object):

    def __init__(self, nRings, pitch):
        self.__nRings = nRings
        self.__nbHexagonalCells = 3 * self.__nRings * (self.__nRings-1) + 1
        self.__pitch = pitch
        self.__hexagonalIndicesPerRow = [None]*(2*self.__nRings-1)
        self.__computeHexagonalIndicesInRows()
        self.__connectivityDict = {}
        self.__computeConnectivityFromHexagonalIndexInRows()
        # Theoretical connectivies =
        # cells away from the boundary layer have 6 connectivities
        # each. The number of such cells correspond to the number of
        # cells in a HoneyCombMesh with nRings - 1.
        # Cells on the boundary but not on the tips have 4
        # connectivities each. There are 6 times (each side) nRings -
        # 2 such cells.
        # Cells on the tips have 3 connectivities each and there are 6
        # of them
        nAuxMeshCells = 3 * (self.__nRings-1) * (self.__nRings-1 - 1) + 1
        theoNoOfConnex = nAuxMeshCells * 6 + 6 * (nRings - 2) * 4 + 6 * 3
        print("Theoretical number of connex = ", theoNoOfConnex)
        print("Number of connex in map = ", len(self.__connectivityDict.keys()))
        assert(theoNoOfConnex == len(self.__connectivityDict.keys()))

        self.__cellIndicesLists = self.__buildCellIndicesListsForAllSextants()
        self.__cellCentres = self.__computeCellCentres()

    @property
    def pitch(self):
        return self.__pitch

    @property
    def halfPitch(self):
        return 0.5 * self.__pitch

    @property
    def halfSide(self):
        return 0.5 * self.__pitch / (3**0.5)
    
    @property
    def domain(self):
        self.__domain

    def getHexagonalIndicesPerRow(self):
        return self.__hexagonalIndicesPerRow

    def __setHexagonalIndicesPerRow(self, hexagonalIndicesPerRow):
        self.__hexagonalIndicesPerRow = hexagonalIndicesPerRow

    def getConnectivityMap(self):
        return self.__connectivityDict

    def __setConnectivityMap(self, connectivityDict):
        self.__connectivityDict = connectivityDict

    def getNumberOfCells(self):
        return self.__nbHexagonalCells

    def __setNumberOfCells(self, nbHexagonalCells):
        self.__nbHexagonalCells = nbHexagonalCells

    def getCellIndicesLists(self):
        return self.__cellIndicesLists

    def __setCellIndicesLists(self, cellIndicesLists):
        self.__cellIndicesLists = cellIndicesLists

    def getCellCentres(self):
        cellCentresDict = {}
        for i in range(self.__nbHexagonalCells):
            cellCentresDict[i] = self.__cellCentres[i]
        return cellCentresDict

    def __setCellCentres(self, cellCentres):
        self.__cellCentres = cellCentres

    def __computeHexagonalIndicesInRows(self):
        nbHexPerCol = 0
        nbBottomCol = self.__nRings
        nbTopCol = 0
        currentHexIndex = 0
        
        for rowIndex in range(self.__nRings):
            # Count from bottom layer to middle layer
            indexInRow = [None] * (nbBottomCol + nbTopCol)
            nbHexPerCol = self.__nRings
            currentHexIndex = self.__nRings - rowIndex - 1
            for colIndex in range(nbBottomCol):
                indexInRow[nbBottomCol + nbTopCol - colIndex - 1] = currentHexIndex
                nbHexPerCol += 1
                currentHexIndex += nbHexPerCol
            currentHexIndex -= nbHexPerCol
            nbHexPerCol -= 1
            currentHexIndex += nbHexPerCol
            for colIndex in range(nbTopCol):
                indexInRow[nbTopCol - colIndex - 1] = currentHexIndex
                nbHexPerCol -= 1
                currentHexIndex += nbHexPerCol
            self.__hexagonalIndicesPerRow[rowIndex] = indexInRow
            nbTopCol += 1

        nbTopCol = self.__nRings
        nbBottomCol -= 1
        initialHexIndex = currentHexIndex - nbHexPerCol

        for rowIndex in range(self.__nRings - 1):
            nbBottomCol -= 1
            nbHexPerCol = self.__nRings
            initialHexIndex -= 1
            currentHexIndex = initialHexIndex
            indexInRow = [None] * (nbTopCol + nbBottomCol)
            for colIndex in range(nbTopCol):
                indexInRow[colIndex] = currentHexIndex
                nbHexPerCol += 1
                currentHexIndex -= nbHexPerCol
            currentHexIndex += nbHexPerCol
            nbHexPerCol -= 1
            currentHexIndex -= nbHexPerCol
            for colIndex in range(nbBottomCol):
                indexInRow[nbTopCol+colIndex] = currentHexIndex
                nbHexPerCol -= 1
                currentHexIndex -= nbHexPerCol
            self.__hexagonalIndicesPerRow[self.__nRings + rowIndex] = indexInRow

    def __computeConnectivityFromHexagonalIndexInRows(self):
        if self.__nRings > 1:
            for rowIndex in range(self.__nRings):
                currentRow = self.__hexagonalIndicesPerRow[rowIndex]
                for colIndex in range(len(currentRow)):
                    currentHexIndex = currentRow[colIndex]
                    if rowIndex > 0:  # To consider only hexagons away from boundary layer
                        bottomRow = self.__hexagonalIndicesPerRow[rowIndex - 1]
                        if colIndex < len(bottomRow):
                            # for current haxagon and bottom edge -> associate neighbour from bottom layer and its edge
                            currentTuple = (currentHexIndex, 0)
                            neighbourTuple = (bottomRow[colIndex], 3)
                            # for current haxagon and SE edge -> associate neighbour from bottom layer and its edge
                            self.__connectivityDict[currentTuple] = neighbourTuple
                        if colIndex > 0:
                            currentTuple = (currentHexIndex, 1)
                            neighbourTuple = (bottomRow[colIndex - 1], 4)
                            self.__connectivityDict[currentTuple] = neighbourTuple
                    if colIndex > 0:
                        currentTuple = (currentHexIndex, 2)
                        neighbourTuple = (currentRow[colIndex-1], 5)
                        self.__connectivityDict[currentTuple] = neighbourTuple
                    if rowIndex < len(self.__hexagonalIndicesPerRow[rowIndex]) - 1:
                        topRow = self.__hexagonalIndicesPerRow[rowIndex + 1]
                        if colIndex < len(topRow):
                            currentTuple = (currentHexIndex, 3)
                            if rowIndex == self.__nRings - 1 and colIndex > 0:
                                # at this frontier, topRow is 1 cell shorter than currentRow
                                neighbourTuple = (topRow[colIndex-1], 0)
                                self.__connectivityDict[currentTuple] = neighbourTuple
                            elif rowIndex == self.__nRings - 1 and colIndex == 0:  # to deal with (last cell, 3)
                                pass
                            else:
                                # at this frontier, topRow is 1 cell longer than currentRow
                                neighbourTuple = (topRow[colIndex], 0)
                                self.__connectivityDict[currentTuple] = neighbourTuple
                        if colIndex == len(topRow):  # to deal with (cell 0, 3)
                            currentTuple = (currentHexIndex, 3)
                            neighbourTuple = (topRow[colIndex-1], 0)
                            self.__connectivityDict[currentTuple] = neighbourTuple
                        if colIndex < len(topRow) - 1:
                            currentTuple = (currentHexIndex, 4)
                            if rowIndex == self.__nRings - 1:
                                # at this frontier, topRow is 1 cell shorter than currentRow
                                neighbourTuple = (topRow[colIndex], 1)
                            else:
                                # at this frontier, topRow is 1 cell longer than currentRow
                                neighbourTuple = (topRow[colIndex + 1], 1)
                            self.__connectivityDict[currentTuple] = neighbourTuple
                        if colIndex == len(topRow) - 1:
                            currentTuple = (currentHexIndex, 4)
                            neighbourTuple = (topRow[colIndex], 1)
                            self.__connectivityDict[currentTuple] = neighbourTuple
                    if colIndex < len(currentRow) - 1:
                        currentTuple = (currentHexIndex, 5)
                        neighbourTuple = (currentRow[colIndex + 1], 2)
                        self.__connectivityDict[currentTuple] = neighbourTuple

            for rowIndex in range(self.__nRings, 2 * self.__nRings - 1):
                currentRow = self.__hexagonalIndicesPerRow[rowIndex]
                for colIndex in range(len(currentRow)):
                    currentHexIndex = currentRow[colIndex]
                    if rowIndex > 0:  # To consider only hexagons away from boundary layer
                        bottomRow = self.__hexagonalIndicesPerRow[rowIndex - 1]
                        # for current haxagon and bottom edge -> associate neighbour from bottom layer and its edge
                        currentTuple = (currentHexIndex, 0)
                        neighbourTuple = (bottomRow[colIndex + 1], 3)
                        self.__connectivityDict[currentTuple] = neighbourTuple
                        currentTuple = (currentHexIndex, 1)
                        neighbourTuple = (bottomRow[colIndex], 4)
                        self.__connectivityDict[currentTuple] = neighbourTuple
                    if colIndex > 0:
                        currentTuple = (currentHexIndex, 2)
                        neighbourTuple = (currentRow[colIndex - 1], 5)
                        self.__connectivityDict[currentTuple] = neighbourTuple
                    if rowIndex < len(self.__hexagonalIndicesPerRow)-1:
                        topRow = self.__hexagonalIndicesPerRow[rowIndex + 1]
                        if colIndex > 0:  # to neglect all cells on the top right border of the geometry
                            currentTuple = (currentHexIndex, 3)
                            neighbourTuple = (topRow[colIndex-1], 0)
                            self.__connectivityDict[currentTuple] = neighbourTuple
                        if colIndex < len(topRow) - 1:
                            currentTuple = (currentHexIndex, 4)
                            neighbourTuple = (topRow[colIndex], 1)
                            self.__connectivityDict[currentTuple] = neighbourTuple
                        if colIndex == len(topRow) - 1:
                            currentTuple = (currentHexIndex, 4)
                            neighbourTuple = (topRow[colIndex], 1)
                            self.__connectivityDict[currentTuple] = neighbourTuple
                    if colIndex < len(currentRow) - 1:
                        currentTuple = (currentHexIndex, 5)
                        neighbourTuple = (currentRow[colIndex + 1], 2)
                        self.__connectivityDict[currentTuple] = neighbourTuple
        else:
            warnings.warn("Connectivity map is empty since mesh has only one hexagon", RuntimeWarning)

    def __buildCellIndicesListsForAllSextants(self):
        if self.__nbHexagonalCells == 1:
            return [[0]] * 6

        cellIndicesLists = [None]*6  # for all sextants

        # first sextant
        currentListForSextant = []
        for hexagonIndex in range(self.__nbHexagonalCells):
            currentListForSextant.append(hexagonIndex)
        cellIndicesLists[0] = currentListForSextant

        # second sextant
        currentListForSextant = []

        nbHexagonsPerLine = self.__nRings
        firstHexagonInLineIndex = self.__nRings - 1
        for lineIndex in range(-self.__nRings, 0):
            for hexagonIndex in range(nbHexagonsPerLine):
                currentListForSextant.append(firstHexagonInLineIndex - hexagonIndex)
            nbHexagonsPerLine += 1
            firstHexagonInLineIndex += nbHexagonsPerLine

        firstHexagonInLineIndex -= nbHexagonsPerLine
        nbHexagonsPerLine -= 2
        firstHexagonInLineIndex += nbHexagonsPerLine
        for lineIndex in range(1, self.__nRings):
            for hexagonIndex in range(nbHexagonsPerLine):
                currentListForSextant.append(firstHexagonInLineIndex - hexagonIndex)
            nbHexagonsPerLine -= 1
            firstHexagonInLineIndex += nbHexagonsPerLine
        cellIndicesLists[1] = currentListForSextant

        # third sextant
        currentListForSextant = []
        for rowIndex in range(2 * self.__nRings - 1):
            currentRow = self.__hexagonalIndicesPerRow[rowIndex]
            for colIndex in range(len(currentRow)):
                hexagonIndex = currentRow[colIndex]
                currentListForSextant.append(hexagonIndex)
        cellIndicesLists[2] = currentListForSextant

        # fourth, fifth and sixth sextant
        cellIndicesLists[3] = [None] * self.__nbHexagonalCells
        cellIndicesLists[4] = [None] * self.__nbHexagonalCells
        cellIndicesLists[5] = [None] * self.__nbHexagonalCells
        currentCellIndex = 0
        for i in range(self.__nbHexagonalCells - 1, -1, -1):
            cellIndicesLists[3][currentCellIndex] = cellIndicesLists[0][i]
            cellIndicesLists[4][currentCellIndex] = cellIndicesLists[1][i]
            cellIndicesLists[5][currentCellIndex] = cellIndicesLists[2][i]
            currentCellIndex += 1

        return cellIndicesLists

    def __computeCellCentres(self):
        # cell centres are appended in the order of the number from 0,
        # 1, ..., nHex
        cellCentres = []

        # u1, u2 and u3 are the directional vectors following the
        # directions of the hexagon
        u1 = (0.5*self.__pitch*3**0.5, -0.5*self.__pitch)
        u2 = (0.5*self.__pitch*3**0.5, 0.5*self.__pitch)
        u3 = (self.__pitch*0., self.__pitch*1.)

        # loop for cells below median line inclusive
        for i in range(-self.__nRings+1, 1):
            x0 = u3[0]*i - u1[0] * (self.__nRings - 1)
            y0 = u3[1]*i - u1[1] * (self.__nRings - 1)
            numberOfCellsPerLine = 2 * self.__nRings - 1 + i
            for j in range(numberOfCellsPerLine):
                cellCentres.append((x0, y0))
                x0 += u1[0]
                y0 += u1[1]

        # loop for cells strictly above median line
        for i in range(1, self.__nRings):
            x0 = u2[0]*i - u1[0] * (self.__nRings - 1)
            y0 = u2[1]*i - u1[1] * (self.__nRings - 1)
            numberOfCellsPerLine = 2 * self.__nRings - 1 - i
            for j in range(numberOfCellsPerLine):
                cellCentres.append((x0, y0))
                x0 += u1[0]
                y0 += u1[1]

        return cellCentres

    def buildCellNodes(self, coords):
        p = 0.5 * self.__pitch
        q = 0.5 * self.__pitch / (3**0.5)
        x, y = coords
        nodes = []
        nodes.append( (x+q,   y+p) )
        nodes.append( (x-q,   y+p) )
        nodes.append( (x-2*q, y  ) )
        nodes.append( (x-q  , y-p) )
        nodes.append( (x+q  , y-p) )
        nodes.append( (x+2*q, y  ) )
        # add 1st node again to close polygon
        nodes.append( (x+q,   y+p) )

        return nodes

    def buildCellNodesFromCentre(self, centres=None):
        points = []

        if centres is None:
            centres = self.__cellCentres

        for coords in centres:
            nodes = self.buildCellNodes(coords)
            points.append(nodes)

        return points

    def plotMesh(self, values=None,
                 pts=None,
                 polyHull=None,
                 colormap=None,
                 label=False,
                 limits=None,
                 fontsizes=[10, 20],
                 fmt=None
                 ):
        import matplotlib.pyplot as plt
        import matplotlib
        from matplotlib.patches import Polygon
        from matplotlib.collections import PatchCollection
        from matplotlib import cm

        if pts is None:
            pts = self.buildCellNodesFromCentre()

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        plt.box(on=None)
        patches = []
        for pt in pts:
            hexagon = Polygon(pt, closed=True, fill=False)
            patches.append(hexagon)

        if values is not None:
            if colormap:
                cmap = colormap
            else:
                cmap = matplotlib.cm.bwr
            # set nan to white
            cmap.set_bad(color='white')
            p = PatchCollection(patches, cmap=cmap, edgecolors="black")
            values = np.ma.array (values, mask=np.isnan(values))
            assert(len(patches) == values.size)
            p.set_array(values)
            if label:
                fmt = fmt if fmt else "{:1.1e}"
                for i in range(len(patches)):
                    pt = pts[i]
                    if values[i] != np.nan:
                        ax.annotate(fmt.format(values[i]), ((pt[0][0]+pt[3][0])/2, (pt[0][1]+pt[3][1])/2), color='black', weight='bold', fontsize=fontsizes[0], ha='center', va='center')
            cbar = plt.colorbar(p)
            for t in cbar.ax.get_yticklabels():
                t.set_fontsize(fontsizes[1])
            if limits:
                assert(len(limits) == 2)
                p.set_clim([limits[0], limits[1]])
        else:
            plt.axis('off')
            p = PatchCollection(patches, edgecolors="black", facecolors="white")
            if colormap:
                colors = 100*np.random.rand(len(patches))
                p.set_array(np.array(colors))
            if label:
                fmt = fmt if fmt else "{:d}"
                for i in range(len(patches)):
                    pt = pts[i]
                    ax.annotate(fmt.format(i), ((pt[0][0]+pt[3][0])/2, (pt[0][1]+pt[3][1])/2), color='black', weight='bold', fontsize=6, ha='center', va='center')
            if polyHull:
                hull = [c for c in polyHull]
                hull.append(polyHull[0])
                xs, ys = zip(*hull)
                plt.plot(xs, ys)

                ghostCells = self.getGhostCells()
                if ghostCells:
                    ghostPts = self.buildCellNodesFromCentre(list(self.getGhostCellCentres().values()))
                    ghostPatches = []
                    for pt in ghostPts:
                        hexagon = Polygon(pt, closed=True, fill=False)
                        ghostPatches.append(hexagon)
                    gp = PatchCollection(ghostPatches,
                                         edgecolors="black",
                                         facecolors="white",
                                         linestyles="dotted"
                                         )
                    if label:
                        fmt = fmt if fmt else "{:d}"
                        nCells = len(patches)
                        for i in range(len(ghostPatches)):
                            pt = ghostPts[i]
                            index = nCells + i
                            ax.annotate(fmt.format(index), ((pt[0][0]+pt[3][0])/2, (pt[0][1]+pt[3][1])/2), color='black', weight='bold', fontsize=6, ha='center', va='center')
                    ax.add_collection(gp)

        ax.add_collection(p)
        ax.axis('scaled')  # ax.autoscale(enable=True)

        plt.tight_layout()
        plt.show()