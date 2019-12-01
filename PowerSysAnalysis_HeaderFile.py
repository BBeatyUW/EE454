"""
Power Flow Analysis: Support Functions
Created By: 
    Brandon Beaty
    Eliot Nichols
"""

import numpy as np
from numpy.linalg import inv
import pandas as pd


""" 
Imports Bus and line data from excel sheets
Takes in an array containing ['File Location', 'Sheet Name']
Returns two panda data frames for the bus and line data
"""
def import_BusAndLineData(BusData_Location, LineData_Location):
    BusData = pd.read_excel(BusData_Location[0], sheet_name=BusData_Location[1])
    LineData = pd.read_excel(LineData_Location[0], sheet_name=LineData_Location[1])
    return BusData, LineData


"""
Builds G and B matrices to be used in Power Flow calculations
Takes in data frame containing all line information, and number of busses in system
Returns G and B arrays
"""
def build_AdmittanceMatrix(LineData, sys_Size):
    col = np.array(LineData.columns)
    line_From = np.array(LineData[col[0]])
    line_To = np.array(LineData[col[1]])
    line_R = np.array(LineData[col[2]])
    line_X = np.array(LineData[col[3]])
    line_Z = np.array(LineData[col[2]]) + 1j*np.array(LineData[col[3]])
    line_Y = 1/line_Z
    line_B = np.array(LineData[col[4]])
    line_Fmax = np.array(LineData[col[5]])
    sys_Y = np.array([[0 for j in range(sys_Size)] for i in range(sys_Size)], dtype = complex)
    sys_G = np.zeros((sys_Size, sys_Size))
    sys_B = np.zeros((sys_Size, sys_Size))
    
    #X_ij
    for i in range(sys_Size): #Row
        for j in range(sys_Size): #Column
            if i==j: # Diagonal, sum of Y(From==i || To==i) + .5B(From==i || To ==i)
                sys_Y[i][j] = np.sum(line_Y[np.array(line_From==i+1) + np.array(line_To==i+1)]) \
                        +.5j*np.sum(line_B[np.array(line_From==i+1) + np.array(line_To==i+1)])
            elif i<j: #Non Diagonal, -Y(From==i && To==j)
                sys_Y[i][j] = -1*line_Y[np.array(line_From==i+1) * np.array(line_To==j+1)]
            else: #i>j =[j][i]
               sys_Y[i][j] = sys_Y[j][i] 
    sys_G = sys_Y.real
    sys_B = sys_Y.imag
    return sys_G, sys_B

"""
Parses intial bus information from data
Takes in Bus Data data frame
Returns sys_:
    LoadP - active power consumed at node
    LoadQ - reactive power consumed at node
    BusType - type of bus<(S)lack, (G)enerator, (D)rain>
    PGen - Active Power produced by each generator node
    VRef - Reference voltages at PV busses
"""
def init_BusData(BusData):
    col = np.array(BusData.columns)
    sys_LoadP = np.array(BusData[col[1]])
    sys_LoadQ = np.array(BusData[col[2]])
    sys_BusType = np.array(BusData[col[3]])
    sys_PGen = np.array(BusData[col[4]])
    sys_VRef = np.array(BusData[col[5]])
    return sys_LoadP, sys_LoadQ, sys_BusType, sys_PGen, sys_VRef

"""
Initializes System Data for processing
Takes in sys_:
    LoadP - active power consumed at node
    LoadQ - reactive power consumed at node
    BusType - type of bus<(S)lack, (G)enerator, (D)rain>
    PGen - Active Power produced by each generator node
    VRef - Reference voltages at PV busses
Returns a 2D array containing each node's current information
[0] - Voltage (V)
[1] - Angle (T)
[2] - Active Power (P)
[3] - P(T,V)-P_inj
[4] - Reactive Power (Q)
[5] - Q(T,V)-Q_inj
"""
def init_SysData(sys_LoadP, sys_LoadQ, sys_BusType, sys_PGen, sys_VRef):
    n= sys_LoadP.size
    sys_Data = np.zeros((n,6))
    sys_Data[:,0] = sys_VRef
    sys_Data[:,1] = np.zeros(n)
    sys_Data[:,2] = sys_PGen-sys_LoadP
    sys_Data[:,3] = sys_PGen-sys_LoadP
    sys_Data[:,4] = -sys_LoadQ
    sys_Data[:,5] = -sys_LoadQ
       
    return sys_Data

"""
Determines Jacobian value for a given cell
Takes in: i, j, n, V_i, V_j, T_i, T_j, P_i, Q_i, G_ij, B_ij
Returns Jacobian cell value
"""
def Jacobian_PowerFlow(i, j, n, V_i, V_j, T_i, T_j, P_i, Q_i, G_ij, B_ij):
    if (i<n and j<n ): #J_11 dP/dT
        if(i==j):
            J = -Q_i-B_ij*V_i**2
        else:
            J =  V_i*V_j*(G_ij*np.sin(T_i-T_j)-B_ij*np.cos(T_i-T_j))
    elif(i<n and j>=n ): #J_12 dP/dV
        if(i==j):
            J = (P_i/V_i) + G_ij*V_i
        else:
            J = V_i*(G_ij*np.cos(T_i-T_j)+B_ij*np.sin(T_i-T_j))
    elif(i>=n and j<n ): #J_21 dQ/dT
        if(i==j):
            J = P_i-G_ij*V_i**2
        else:
            J = -V_i*V_j*(G_ij*np.cos(T_i-T_j)+B_ij*np.sin(T_i-T_j))
    else: #(i>=n and j>= n) J_22 dQ/dV
        if(i==j):
            J = (Q_i/V_i)-B_ij*V_i
        else:
            J = V_i*(G_ij*np.sin(T_i-T_j)-B_ij*np.cos(T_i-T_j))
    return J

"""
Processes 1 iteration of current system data
Takes in sys_Data, a 2D array containing each node's current information
[0] - Voltage (V)
[1] - Angle (T)
[2] - Active Power (P) PGen-PLoad
[3] - P(T,V)-P_inj
[4] - Reactive Power (Q) QGen-QLoad
[5] - Q(T,V)-Q_inj
As well as, the systems G and B matrices, and node types
Returns the updated array
"""
def update_SysData(sys_Data, sys_G, sys_B, sys_BusType):
    nodes = sys_BusType.size
    n = nodes-1
    """ Determine Jacobian """
    J = np.zeros((2*n,2*n))
    for i in range(2*n):
        for j in range(2*n):
            J[i,j] = Jacobian_PowerFlow(i, j, 2*n, \
             sys_Data[(i%n)+1 ,0], sys_Data[(j%n)+1 ,0],\
             sys_Data[(i%n)+1 ,1], sys_Data[(j%n)+1 ,1],\
             sys_Data[(i%n)+1 ,3] , sys_Data[(i%n)+1 ,5],\
             sys_G[(i%n)+1, (j%n)+1], sys_B[(i%n)+1, (j%n)+1])
    
    """ Determine inverse of Jacobian """
    J_inv = inv(J)
    
    """ Calculate Delta V's, Theta's, and update """
    Delta = J_inv @ np.append(sys_Data[1:nodes,1], sys_Data[1:nodes,0])
    Delta_T = np.append([0], Delta[0:n])
    Delta_V = np.append([0], Delta[n:2*n])
    sys_Data[:,0] += Delta_V
    sys_Data[:,1] += Delta_T
    
    """ Update Q_inj[4], P(T,V)[3],Q(T,V)[5] """ 
    """Explicit Update: Genator Q_inj, Slack P_inj"""
    for i in range(nodes):
        for j in range(nodes):
            if i==0:
                sys_Data[i,2] += sys_Data[i,0]*sys_Data[j,0]*\
                            (sys_G[i,j]*np.cos(sys_Data[i,1]-sys_Data[j,1])+\
                             sys_B[i,j]*np.sin(sys_Data[i,1]-sys_Data[j,1]))
            if sys_BusType[i]!='D':
                sys_Data[i,4] += sys_Data[i,0]*sys_Data[j,0]*\
                            (sys_G[i,j]*np.sin(sys_Data[i,1]-sys_Data[j,1])-\
                             sys_B[i,j]*np.cos(sys_Data[i,1]-sys_Data[j,1]))
                            
    """Implicit Update: P(T,V)-P_inj,Q(T,V)-Q_inj"""
    for i in range(nodes):
        sys_Data[i,3] = -1*sys_Data[i,2]
        sys_Data[i,5] = -1*sys_Data[i,4]
        for j in range(nodes):
            sys_Data[i,3] += sys_Data[i,0]*sys_Data[j,0]*\
                            (sys_G[i,j]*np.cos(sys_Data[i,1]-sys_Data[j,1])+\
                             sys_B[i,j]*np.sin(sys_Data[i,1]-sys_Data[j,1]))
            sys_Data[i,5] += sys_Data[i,0]*sys_Data[j,0]*\
                        (sys_G[i,j]*np.sin(sys_Data[i,1]-sys_Data[j,1])-\
                         sys_B[i,j]*np.cos(sys_Data[i,1]-sys_Data[j,1]))
                        
    return sys_Data


"""
Code for testing purposes
"""


"""
print(sys_G)
print()
print(sys_B)
"""

"""
dict_BusData = {df_BusData.columns[0]:df_BusData[df_BusData.columns[0]]}

print(df_BusData['P MW'][df_BusData['Type']=='G'])

line_From = np.array(df_LineData['From'])
line_To = np.array(df_LineData['To'])
busNum = np.array(df_BusData['Bus #'])      
ind = np.array(df_LineData.columns)
line_R = np.array(df_LineData[ind[2]])
"""

"""
ind1 = np.array(line_From==1)
ind2 = np.array(line_To==2)
ind3 = ind1 * ind2
print(-1*line_R[np.array(line_From==1) * np.array(line_To==2)])

for i in range(df_BusData.shape[0]):
    print(i)
#print(df_BusData)
#print(ind.size)
#print(build_AdmittanceMatrix(df_LineData))
"""
"""
test = np.array([[2*i+j for j in range(2)] for i in range(2)], dtype = complex)
print(test)
print(test[:,0].size)
test2 = np.zeros((2,3))
test2[:,1] = np.array([4,5])
print(test2)
#"""

"""
test1 = 0.5*np.ones((2,2))
test2 = 1/test1
#print(test2)
"""
"""
testshape = np.array([[1,2,3],[1,2,3]])
print(testshape.shape)
"""