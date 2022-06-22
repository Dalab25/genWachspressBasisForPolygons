import numpy as np
from honeycomb import HoneyCombMesh
from irregular_pentagon import IrregularPentagon
from geometrical_elements import hexCycle, pentaCycle
from hexagonal_functions_generalised_basis import HexagonalFunctionsGeneralisedBasis
from pentagonal_functions_generalised_basis import PentagonalFunctionsGeneralisedBasis

def getWSF(p, order, nSides, boolHomothetic=True, bool_gout=False):
    """
    Get the WSF functions for a specific order and a given polygon
    """
    # pitch = distance between two parallel sides
    # p = 1.

    # create hexagonal mesh with 1 hexagonal cell
    if(nSides == 5):
        mesh = IrregularPentagon()
    elif(nSides == 6):
        mesh = HoneyCombMesh(1, p)
    else:
        print("Mesh does not exist for this number of sides")

    # order = 3
    # Build the WSF object to generate shape functions
    if(boolHomothetic):
        coeffsFile = "../results/polygon_" + str(nSides) + "sides/order" + str(order) + "/symbolic_coefficients/solution_order_" + str(order) + "_pitch_1_free_variables_0.txt"
        if(nSides == 5):
            wsf1 = PentagonalFunctionsGeneralisedBasis(order, mesh, p, coeffsFile,homothetic=True)
        
        elif(nSides == 6):
            wsf1 = HexagonalFunctionsGeneralisedBasis(mesh, order, coeffsFile,homothetic=True)
        else:
            print("WSF basis does not exist for this polygon")
    else:
        coeffsFile = "../results/polygon_6sides/order" + str(order) + "/symbolic_coefficients/solution_order_" + str(order) + "_pitch_" + str(p) + "_free_variables_0.txt"
        wsf1 = HexagonalFunctionsGeneralisedBasis(mesh, order, coeffsFile,homothetic=False)
   
    # Get coordinates and geometrical objects
    coordsDict, linesDict = wsf1.getCoordsLines()

    # Choose the appropriate method to build the shape functions for a
    # given order
    if order == 1:
        func = wsf1.computeFirstOrderFunction
    elif order == 2:
        #func = wsf1.computeSecondOrderFunctionGout
        func = wsf1.computeSecondOrderFunction
    elif order == 3:
        if(bool_gout == True):
            func = wsf1.computeThirdOrderFunctionGout
            print("Using Gout basis")
        else:
            func = wsf1.computeThirdOrderFunction
            print("Using basis generated from the symbolic computation")
    elif order == 4:
        func = wsf1.computeFourthOrderFunction
    elif order == 5:
        func = wsf1.computeFifthOrderFunction
    elif order == 6:
        func = wsf1.computeSixthOrderFunction
    else:
        raise ValueError("Order %s is not yet available" % order)

    return func, coordsDict, linesDict

def constraints_verification(order, p, nSides, boolHomothetic, bool_gout=False):
    """
    Verification of wi(a_j) = delta_ij
    """
    func, coords, lines = getWSF(p, order, nSides, boolHomothetic, bool_gout)
    
    L_w_i = []
    L_w_i_1 = []
    L_w_i_2 = []
    L_w_i_3 = []
    L_w_i_4 = []
    
    funcDict = {}

    if(order == 1):
        
        for i in range(1,7):
            
            wi = func(i)
            funcDict['w%s' % i] = wi
            ajx, ajy= coords['a%d' % i][0], coords['a%d' %i][1]
            
            err = 1 - wi(ajx, ajy)
            L_w_i.append(err)
            
    elif(order == 2):
        for i in range(1,7):

            #getting the wi, wi_j functions
            wi,  wiip1= func(i)
            
            funcDict['w%s' % i] = wi
            nextI = hexCycle(i+1)
            funcDict['w%s%s' % (i,nextI)] = wiip1
            
            
            ajx, ajy= coords['a%d' % i][0], coords['a%d' %i][1]
            ajjp1x, ajjp1y = coords['a%d%d' % (i,nextI)][0], coords['a%d%d' % (i,nextI)][1]
            
            if(i==1):
                aj_1jx, aj_1jy = coords['a%d%d' % (6,i)][0], coords['a%d%d' % (6,i)][1]
            else:    
                aj_1jx, aj_1jy = coords['a%d%d' % (i-1,i)][0], coords['a%d%d' % (i-1,i)][1]
                
            err = 1 - wi(ajx, ajy)
            L_w_i.append(err)
            L_w_i.append(wi(ajjp1x, ajjp1y))
            L_w_i.append(wi(aj_1jx, aj_1jy))  
            
            err = 1 - wiip1(ajjp1x, ajjp1y)
            L_w_i_1.append(err)
            
    
    elif order == 3:
        
        for i in range(1,nSides + 1):
            if(nSides == 5):
                i -= 1
            #getting the wi, wi_j functions
            wi,  wiiip1, wip1ip1i= func(i)
            
            funcDict['w%s' % i] = wi
            if(nSides == 5):
                nextI = pentaCycle(i + 1)
            else:
                nextI = hexCycle(i + 1)
            funcDict['w%s%s%s' % (nextI, nextI, i)] = wip1ip1i
            funcDict['w%s%s%s' % (i, i, nextI)] = wiiip1
            
            
            ajx, ajy= coords['a%d' % i][0], coords['a%d' %i][1]
            ajjjp1x, ajjjp1y = coords['a%d%d%d' % (i,i,nextI)][0], coords['a%d%d%d' % (i,i,nextI)][1]
            ajp1jp1jx, ajp1jp1jy = coords['a%d%d%d' % (nextI,nextI,i)][0], coords['a%d%d%d' % (nextI,nextI,i)][1]
            
            if(nSides == 5):
                if(i==0):
                    aj_1j_1jx, aj_1j_1jy = coords['a%d%d%d' % (4,4,i)][0], coords['a%d%d%d' % (4,4,i)][1]
                    ajjj_1x, ajjj_1y = coords['a%d%d%d' % (i,i,4)][0], coords['a%d%d%d' % (i,i,4)][1]
                else:    
                    aj_1j_1jx, aj_1j_1jy = coords['a%d%d%d' % (i-1,i-1,i)][0], coords['a%d%d%d' % (i-1,i-1,i)][1]
                    ajjj_1x, ajjj_1y = coords['a%d%d%d' % (i,i,i-1)][0], coords['a%d%d%d' % (i,i,i-1)][1]
               
            else:
                if(i==1):
                    aj_1j_1jx, aj_1j_1jy = coords['a%d%d%d' % (nSides,nSides,i)][0], coords['a%d%d%d' % (nSides,nSides,i)][1]
                    ajjj_1x, ajjj_1y = coords['a%d%d%d' % (i,i,nSides)][0], coords['a%d%d%d' % (i,i,nSides)][1]
                else:    
                    aj_1j_1jx, aj_1j_1jy = coords['a%d%d%d' % (i-1,i-1,i)][0], coords['a%d%d%d' % (i-1,i-1,i)][1]
                    ajjj_1x, ajjj_1y = coords['a%d%d%d' % (i,i,i-1)][0], coords['a%d%d%d' % (i,i,i-1)][1]
                    
            err = 1 - wi(ajx, ajy)
            L_w_i.append(err)
            L_w_i.append(wi(ajjjp1x, ajjjp1y))
            L_w_i.append(wi(ajp1jp1jx, ajp1jp1jy))
            L_w_i.append(wi(aj_1j_1jx, aj_1j_1jy))  
            L_w_i.append(wi(ajjj_1x, ajjj_1y))
            
            err = 1 - wiiip1(ajjjp1x, ajjjp1y)
            L_w_i_1.append(err)
            L_w_i_1.append(wiiip1(ajp1jp1jx, ajp1jp1jy)) 

            err = 1 - wip1ip1i(ajp1jp1jx, ajp1jp1jy)
            L_w_i_2.append(err)
            L_w_i_2.append(wip1ip1i(ajjjp1x, ajjjp1y))   
            
       
    elif order == 4:
        
        for j in range(1, nSides + 1):
            if(nSides == 5):
                j -= 1
            wi,wiiip1,wiip1,wip1ip1i = func(j)
            
            funcDict['w%s' % j] = wi
            if(nSides == 5):
                nextJ = pentaCycle(j + 1)
            else:
                nextJ = hexCycle(j + 1)
                
            funcDict['w%s%s%s'%(j,j,nextJ)] = wiiip1
            funcDict['w%s%s'% (j,nextJ)] = wiip1
            funcDict['w%s%s%s'%(nextJ,nextJ,j)] = wip1ip1i
            
            ajx, ajy= coords['a%d' % j][0], coords['a%d' %j][1]
            ajjjp1x, ajjjp1y = coords['a%d%d%d' % (j,j,nextJ)][0], coords['a%d%d%d' % (j,j,nextJ)][1]
            ajjp1x, ajjp1y = coords['a%d%d' % (j,nextJ)][0], coords['a%d%d' % (j,nextJ)][1]
            ajp1jp1jx, ajp1jp1jy = coords['a%d%d%d' % (nextJ,nextJ,j)][0], coords['a%d%d%d' % (nextJ,nextJ,j)][1]
                
            if(nSides == 5):
                if(j == 0):
                    aj_1j_1jx, aj_1j_1jy = coords['a%d%d%d' % (4,4,j)][0], coords['a%d%d%d' % (4,4,j)][1]
                    aj_1jx, aj_1jy = coords['a%d%d' % (4,j)][0], coords['a%d%d' % (4,j)][1]
                    ajjj_1x, ajjj_1y = coords['a%d%d%d' % (j,j,4)][0], coords['a%d%d%d' % (j,j,4)][1]
                else:    
                    aj_1j_1jx, aj_1j_1jy = coords['a%d%d%d' % (j-1,j-1,j)][0], coords['a%d%d%d' % (j-1,j-1,j)][1]
                    aj_1jx, aj_1jy = coords['a%d%d' % (j-1,j)][0], coords['a%d%d' % (j-1,j)][1]
                    ajjj_1x, ajjj_1y = coords['a%d%d%d' % (j,j,j-1)][0], coords['a%d%d%d' % (j,j,j-1)][1]
            else:
                if(j==1):
                    aj_1j_1jx, aj_1j_1jy = coords['a%d%d%d' % (6,6,j)][0], coords['a%d%d%d' % (6,6,j)][1]
                    aj_1jx, aj_1jy = coords['a%d%d' % (6,j)][0], coords['a%d%d' % (6,j)][1]
                    ajjj_1x, ajjj_1y = coords['a%d%d%d' % (j,j,6)][0], coords['a%d%d%d' % (j,j,6)][1]
                else:    
                    aj_1j_1jx, aj_1j_1jy = coords['a%d%d%d' % (j-1,j-1,j)][0], coords['a%d%d%d' % (j-1,j-1,j)][1]
                    aj_1jx, aj_1jy = coords['a%d%d' % (j-1,j)][0], coords['a%d%d' % (j-1,j)][1]
                    ajjj_1x, ajjj_1y = coords['a%d%d%d' % (j,j,j-1)][0], coords['a%d%d%d' % (j,j,j-1)][1]
                     
            err = 1 - wi(ajx, ajy)
            L_w_i.append(err)            
            L_w_i.append(wi(ajjjp1x, ajjjp1y)) 
            L_w_i.append(wi(ajjp1x, ajjp1y)) 
            L_w_i.append(wi(ajp1jp1jx, ajp1jp1jy))  
            L_w_i.append(wi(aj_1j_1jx, aj_1j_1jy)) 
            L_w_i.append(wi(aj_1jx, aj_1jy))    
            L_w_i.append(wi(aj_1jx, aj_1jy))    

            err = 1 - wiiip1(ajjjp1x, ajjjp1y)
            L_w_i_1.append(err)
            L_w_i_1.append(wiiip1(ajjp1x, ajjp1y))     
            L_w_i_1.append(wiiip1(ajp1jp1jx, ajp1jp1jy))

            err = 1 - wiip1(ajjp1x, ajjp1y)
            L_w_i_2.append(err)
            L_w_i_2.append(wiip1(ajjjp1x, ajjjp1y))     
            L_w_i_2.append(wiip1(ajp1jp1jx, ajp1jp1jy))
            
            err = 1 - wip1ip1i(ajp1jp1jx, ajp1jp1jy)
            L_w_i_3.append(err)
            L_w_i_3.append(wip1ip1i(ajjjp1x, ajjjp1y))
            L_w_i_3.append(wip1ip1i(ajjp1x, ajjp1y))
            
    elif order == 5:            
            
        for j in range(1,7):
            
            wi,wi_1,wi_2,wi_3,wi_4 = func(j)
            funcDict['w%s'% j] = wi
            nextI = hexCycle(j+1)
            funcDict['w%s_1'%j] = wi_1
            funcDict['w%s_2'%j] = wi_2
            funcDict['w%s_3'%j] = wi_3
            funcDict['w%s_4'%j] = wi_4
            
            
            ajx, ajy= coords['a%d' % j][0], coords['a%d' %j][1]
            aj1_x, aj1_y = coords['a%d1' %j][0], coords['a%d1' %j][1]
            aj2_x, aj2_y = coords['a%d2' %j][0], coords['a%d2' %j][1]
            aj3_x, aj3_y = coords['a%d3' %j][0], coords['a%d3' %j][1]
            aj4_x, aj4_y = coords['a%d4' %j][0], coords['a%d4' %j][1]

            if(j==1):
                aj_1_1_x, aj_1_1_y = coords['a61'][0], coords['a61'][1]
                aj_1_2_x, aj_1_2_y = coords['a62'][0], coords['a62'][1]
                aj_1_3_x, aj_1_3_y = coords['a63'][0], coords['a63'][1]
                aj_1_4_x, aj_1_4_y = coords['a64'][0], coords['a64'][1]
                
            else:    
                j_1 = j - 1
                aj_1_1_x, aj_1_1_y = coords['a%d1' %j_1][0], coords['a%d1' %j_1][1]
                aj_1_2_x, aj_1_2_y = coords['a%d1' %j_1][0], coords['a%d1' %j_1][1]
                aj_1_3_x, aj_1_3_y = coords['a%d1' %j_1][0], coords['a%d1' %j_1][1]
                aj_1_4_x, aj_1_4_y = coords['a%d1' %j_1][0], coords['a%d1' %j_1][1]
                
            err = 1 - wi(ajx, ajy)
            L_w_i.append(err)
            L_w_i.append(wi(aj1_x, aj1_y))
            L_w_i.append(wi(aj2_x, aj2_y))
            L_w_i.append(wi(aj3_x, aj3_y))
            L_w_i.append(wi(aj4_x, aj4_y))
            L_w_i.append(wi(aj_1_1_x, aj_1_1_y))
            L_w_i.append(wi(aj_1_2_x, aj_1_2_y))
            L_w_i.append(wi(aj_1_3_x, aj_1_3_y))
            L_w_i.append(wi(aj_1_4_x, aj_1_4_y))
             
            err = 1 - wi_1(aj1_x, aj1_y)
            L_w_i_1.append(err)
            L_w_i_1.append(wi_1(aj2_x, aj2_y))
            L_w_i_1.append(wi_1(aj3_x, aj3_y))
            L_w_i_1.append(wi_1(aj4_x, aj4_y))
            
            err = 1 - wi_2(aj2_x, aj2_y)
            L_w_i_2.append(err)
            L_w_i_2.append(wi_2(aj1_x, aj1_y))
            L_w_i_2.append(wi_2(aj3_x, aj3_y))
            L_w_i_2.append(wi_2(aj4_x, aj4_y))
            
            err = 1 - wi_3(aj3_x, aj3_y)
            L_w_i_3.append(err)
            L_w_i_3.append(wi_3(aj1_x, aj1_y))
            L_w_i_3.append(wi_3(aj2_x, aj2_y))
            L_w_i_3.append(wi_3(aj4_x, aj4_y))
            
            err = 1 - wi_4(aj4_x, aj4_y)
            L_w_i_4.append(err)
            L_w_i_4.append(wi_4(aj1_x, aj1_y))
            L_w_i_4.append(wi_4(aj2_x, aj2_y))
            L_w_i_4.append(wi_4(aj3_x, aj3_y))
        
    else:
        raise ValueError("Order not implemented")
    
            
    print("Erreur sur les wi", np.linalg.norm(L_w_i)) 
    assert(np.linalg.norm(L_w_i) < 10**-10)
    print("Erreur sur les wi_1", np.linalg.norm(L_w_i_1))
    assert(np.linalg.norm(L_w_i_1) < 10**-10)   
    if(order >= 3):
        print("Erreur sur les wi_2", np.linalg.norm(L_w_i_2)) #ce coeff n'apparait pas à l'ordre 4
        assert(np.linalg.norm(L_w_i_2) < 10**-10)
    if(order >= 4):
        print("Erreur sur les wi_3", np.linalg.norm(L_w_i_3))
        assert(np.linalg.norm(L_w_i_3) < 10**-10)
    if(order == 5):
        print("Erreur sur les wi_4", np.linalg.norm(L_w_i_4))
        assert(np.linalg.norm(L_w_i_4) < 10**-10)
    
       
def main():
    orders = [3]
    bool_gout = False
    boolHomothetic = True
    nSides = 6
    #p = 10
    p = 1
    for order in orders:
        constraints_verification(order, p, nSides, boolHomothetic, bool_gout)            
            
if __name__ == "__main__":
    main()            
            
            