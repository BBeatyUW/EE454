"""
Power Flow Analysis: Support Functions
Split of PQ and PV buses
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
                sys_Y[i][j] = -np.sum(line_Y[np.multiply(np.array(line_From==i+1), np.array(line_To==j+1))])
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
    sys_BusNum = np.array(BusData[col[0]])
    sys_LoadP = np.array(BusData[col[1]])
    sys_LoadQ = np.array(BusData[col[2]])
    sys_BusType = np.array(BusData[col[3]])
    sys_PGen = np.array(BusData[col[4]])
    sys_VRef = np.array(BusData[col[5]])
    return sys_BusNum, sys_LoadP, sys_LoadQ, sys_BusType, sys_PGen, sys_VRef

"""
Initializes System Data for processing
Takes in sys_:
    LoadP - active power consumed at node
    LoadQ - reactive power consumed at node
    BusType - type of bus<(S)lack, (G)enerator, (D)rain>
    PGen - Active Power produced by each generator node
    VRef - Reference voltages at PV busses
Returns a 2D array containing each node's current information
[0] - Bus #
[1] - Voltage (V)
[2] - Angle (T)
[3] - Active Power (P_inj)
[4] - P(T,V)-P_inj (mismatch)
[5] - Reactive Power (Q_inj)
[6] - Q(T,V)-Q_inj (mismatch)
"""
def init_SysData(sys_BusNum, sys_LoadP, sys_LoadQ, sys_BusType, sys_PGen, sys_VRef, sys_G, sys_B, S_Base):
    n= sys_LoadP.size
    sys_Data = np.zeros((n,7))
    sys_Data[:,0] = sys_BusNum
    sys_Data[:,1] = sys_VRef #Sets initial voltages to provided reference
    sys_Data[:,2] = np.zeros(n) #Sets initial angles to zero
    sys_Data[:,3] = (sys_PGen-sys_LoadP)/S_Base #Sets initial power inject to Bus generation minus load in per unit
    sys_Data[sys_BusType=='S',3] = (np.sum(sys_LoadP)-np.sum(sys_PGen))/S_Base #Sets initial guess for active power required from slack bus
    sys_Data[:,5] = (-sys_LoadQ)/S_Base #Sets initial power inject to Bus generation minus load in per unit  
    sys_Data[sys_BusType=='S',5] = (-np.sum(sys_LoadQ))/S_Base #Sets initial guess for reactive power required from slack bus
    for i in range(n): #Sets initial mismatch to calculated power from (V,T) minus expected inject
        sys_Data[i,4] = -sys_Data[i,3]
        sys_Data[i,6] = -sys_Data[i,5] 
        for j in range(n):
            sys_Data[i,4] += sys_Data[i,1]*sys_Data[j,1]*\
                            (sys_G[i,j]*np.cos(sys_Data[i,2]-sys_Data[j,2])+\
                             sys_B[i,j]*np.sin(sys_Data[i,2]-sys_Data[j,2]))
            sys_Data[i,6] += sys_Data[i,1]*sys_Data[j,1]*\
                        (sys_G[i,j]*np.sin(sys_Data[i,2]-sys_Data[j,2])-\
                         sys_B[i,j]*np.cos(sys_Data[i,2]-sys_Data[j,2]))
       
    return sys_Data

"""
Determines Jacobian value for a given J_11 cell (dP/dT)
Takes in: i, j, n, V_i, V_j, T_i, T_j, P_i, Q_i, G_ij, B_ij
Returns Jacobian cell value
"""
def Jacobian_PowerFlow_11(i, j, V_i, V_j, T_i, T_j, P_i, Q_i, G_ij, B_ij):
    if(i==j):
        J = -Q_i - B_ij*(V_i**2)
    else:
        J =  V_i*V_j*(G_ij*np.sin(T_i-T_j)-B_ij*np.cos(T_i-T_j))
    return J

"""
Determines Jacobian value for a given J_12 cell (dP/dV)
Takes in: i, j, n, V_i, V_j, T_i, T_j, P_i, Q_i, G_ij, B_ij
Returns Jacobian cell value
"""
def Jacobian_PowerFlow_12(i, j, V_i, V_j, T_i, T_j, P_i, Q_i, G_ij, B_ij):
    if(i==j):
        J = (P_i/V_i) + G_ij*V_i
    else:
        J = V_i*(G_ij*np.cos(T_i-T_j)+B_ij*np.sin(T_i-T_j))
    return J

"""
Determines Jacobian value for a given J_21 cell (dQ/dT)
Takes in: i, j, n, V_i, V_j, T_i, T_j, P_i, Q_i, G_ij, B_ij
Returns Jacobian cell value
"""
def Jacobian_PowerFlow_21(i, j, V_i, V_j, T_i, T_j, P_i, Q_i, G_ij, B_ij):
    if(i==j):
        J = P_i-G_ij*(V_i**2)
    else:
        J = -V_i*V_j*(G_ij*np.cos(T_i-T_j)+B_ij*np.sin(T_i-T_j))
    return J

"""
Determines Jacobian value for a given J_22 cell (dQ/dV)
Takes in: i, j, n, V_i, V_j, T_i, T_j, P_i, Q_i, G_ij, B_ij
Returns Jacobian cell value
"""
def Jacobian_PowerFlow_22(i, j, V_i, V_j, T_i, T_j, P_i, Q_i, G_ij, B_ij):
    if(i==j):
        J = (Q_i/V_i)-B_ij*V_i
    else:
        J = V_i*(G_ij*np.sin(T_i-T_j)-B_ij*np.cos(T_i-T_j))
    return J


"""
Processes 1 iteration of current system data
Takes in sys_Data, a 2D array containing each node's current information
[0] - Bus #
[1] - Voltage (V)
[2] - Angle (T)
[3] - Active Power (P_inj)
[4] - P(T,V)-P_inj (mismatch)
[5] - Reactive Power (Q_inj)
[6] - Q(T,V)-Q_inj (mismatch)
As well as, the systems G and B matrices, and node types
Returns the updated array
"""
def update_SysData(sys_Data, sys_G, sys_B, sys_BusType):
    n = sys_BusType.size
    D_index = sys_BusType=='D'
    G_index = sys_BusType=='G'
    S_index = sys_BusType=='S'
    
    """Determine Jacobian"""
    J = np.zeros((2*n,2*n))
    for i in range(n):
        for j in range(n):                  #(i, j, V_i,           V_j,           T_i,           T_j,           P_i(T,V),                    Q_i(T,V),                    G_ij,       B_ij)
            J[i,j] =    Jacobian_PowerFlow_11(i, j, sys_Data[i,1], sys_Data[j,1], sys_Data[i,2], sys_Data[j,2], sys_Data[i,4]+sys_Data[i,3], sys_Data[i,6]+sys_Data[i,5], sys_G[i,j], sys_B[i,j]) 
            J[i,j+n] =  Jacobian_PowerFlow_12(i, j, sys_Data[i,1], sys_Data[j,1], sys_Data[i,2], sys_Data[j,2], sys_Data[i,4]+sys_Data[i,3], sys_Data[i,6]+sys_Data[i,5], sys_G[i,j], sys_B[i,j])
            J[i+n,j] =  Jacobian_PowerFlow_21(i, j, sys_Data[i,1], sys_Data[j,1], sys_Data[i,2], sys_Data[j,2], sys_Data[i,4]+sys_Data[i,3], sys_Data[i,6]+sys_Data[i,5], sys_G[i,j], sys_B[i,j]) 
            J[i+n,j+n] =Jacobian_PowerFlow_22(i, j, sys_Data[i,1], sys_Data[j,1], sys_Data[i,2], sys_Data[j,2], sys_Data[i,4]+sys_Data[i,3], sys_Data[i,6]+sys_Data[i,5], sys_G[i,j], sys_B[i,j])
    
    """Remove non-implicit values"""
    for i in range(n-1,-1,-1):
        if S_index[i]:
            J=np.delete(J, i+n, 0)
            J=np.delete(J, i+n, 1)
            J=np.delete(J, i, 0)
            J=np.delete(J, i, 1)
        elif G_index[i]:
            J=np.delete(J, i+n, 0)
            J=np.delete(J, i+n, 1)
    print(J)
    """Determine Inverse"""
    J_inv = inv(J)
    
    """Determine Delta T,V"""
    PQ = np.concatenate((sys_Data[np.invert(S_index), 4], sys_Data[D_index, 6]))
    Delta = -J_inv @ PQ
    Delta_T = Delta[0:sum(np.invert(S_index))]
    Delta_V = Delta[sum(np.invert(S_index)):sum(np.invert(S_index))+sum(D_index)]
    """Update T for non-slack buses, and V for PQ buses"""
    Delta_T_index = 0
    Delta_V_index = 0
    for i in range(n):
        if G_index[i]:
            sys_Data[i,2] += Delta_T[Delta_T_index]
            Delta_T_index += 1
        elif D_index[i]:
            sys_Data[i,1] += Delta_V[Delta_V_index]
            Delta_V_index += 1
            sys_Data[i,2] += Delta_T[Delta_T_index]
            Delta_T_index += 1
    
    """Update P_inj for slack bus, and Q_inj for non PQ buses"""
    for i in range(n):
        if S_index[i]:#Update Slack P_inj
            sys_Data[i,3] = 0
        if (S_index[i] or G_index[i]):#Update non PQ Q_inj
            sys_Data[i,5] = 0
        for j in range(n):
            if S_index[i]:#Update Slack
                sys_Data[i,3] += sys_Data[i,1]*sys_Data[j,1]*((sys_G[i,j]*np.cos(sys_Data[i,2]-sys_Data[j,2]))+(sys_B[i,j]*np.sin(sys_Data[i,2]-sys_Data[j,2])))
            if (S_index[i] or G_index[i]):#Update non PQ
                sys_Data[i,5] += sys_Data[i,1]*sys_Data[j,1]*((sys_G[i,j]*np.sin(sys_Data[i,2]-sys_Data[j,2]))-(sys_B[i,j]*np.cos(sys_Data[i,2]-sys_Data[j,2])))
    
    """Update mismatch columns"""
    for i in range(n):
        sys_Data[i,4] = -sys_Data[i,3]
        sys_Data[i,6] = -sys_Data[i,5]
        for j in range(n):
            sys_Data[i,4] += sys_Data[i,1]*sys_Data[j,1]*((sys_G[i,j]*np.cos(sys_Data[i,2]-sys_Data[j,2]))+(sys_B[i,j]*np.sin(sys_Data[i,2]-sys_Data[j,2])))
            sys_Data[i,6] += sys_Data[i,1]*sys_Data[j,1]*((sys_G[i,j]*np.sin(sys_Data[i,2]-sys_Data[j,2]))-(sys_B[i,j]*np.cos(sys_Data[i,2]-sys_Data[j,2])))
    
    return sys_Data

# =============================================================================
# 
# def updateI_SysData(sys_Data, sys_G, sys_B, sys_BusType):
#     """Determine PQ buses"""
#     n = sys_BusType.size
#     nI = sys_BusType[sys_BusType != 'S'].size
#     sys_DataI = np.zeros((nI, sys_Data.shape[1]))
#     
#     index = 0
#     for i in range(n):
#         if sys_BusType[i] != 'S':
#             sys_DataI[index, :] = sys_Data[i, :]
#             index += 1
#     """ Determine Jacobian """
#     J = np.zeros((2*nI,2*nI))
#     
#     for i in range (nI):
#         for j in range (nI):                  #(i, j, V_i,               V_j,             T_i,              T_j,              P_i,                         Q_i,             G_ij,        B_ij)
#             J[i,j] =      Jacobian_PowerFlow_11(i, j, sys_DataI[i ,1], sys_DataI[j ,1], sys_DataI[i ,2], sys_DataI[j ,2], sys_DataI[i, 4]+sys_DataI[i, 3], sys_DataI[i, 6]+sys_DataI[i, 5], sys_G[i+1, j+1], sys_B[i+1, j+1])
#             J[i,j+nI] =   Jacobian_PowerFlow_12(i, j, sys_DataI[i ,1], sys_DataI[j ,1], sys_DataI[i ,2], sys_DataI[j ,2], sys_DataI[i, 4]+sys_DataI[i, 3], sys_DataI[i, 6]+sys_DataI[i, 5], sys_G[i+1, j+1], sys_B[i+1, j+1])
#             J[i+nI,j] =   Jacobian_PowerFlow_21(i, j, sys_DataI[i ,1], sys_DataI[j ,1], sys_DataI[i ,2], sys_DataI[j ,2], sys_DataI[i, 4]+sys_DataI[i, 3], sys_DataI[i, 6]+sys_DataI[i, 5], sys_G[i+1, j+1], sys_B[i+1, j+1])
#             J[i+nI,j+nI]= Jacobian_PowerFlow_22(i, j, sys_DataI[i ,1], sys_DataI[j ,1], sys_DataI[i ,2], sys_DataI[j ,2], sys_DataI[i, 4]+sys_DataI[i, 3], sys_DataI[i, 6]+sys_DataI[i, 5], sys_G[i+1, j+1], sys_B[i+1, j+1])
#     #print(J)
#     """ Determine inverse of Jacobian """
#     J_inv = inv(J)
#     
#     """ Calculate Delta V's, Theta's, and update V and Theta """
#     PQ_TV = np.concatenate((sys_DataI[:,4], sys_DataI[:,6]), axis=None)
#     Delta = -J_inv @ PQ_TV
#     Delta_T = Delta[0:nI]
#     Delta_V = Delta[nI:2*nI]
#     for i in range(nI):
#         if (sys_BusType[1:n])[i]=='D':
#             sys_DataI[i,1] += Delta_V[i]
#     sys_DataI[:,2] += Delta_T
#     #print(Delta)
#     
#     """Merge PQ update with rest of system"""
#     index=0
#     for i in range(n):
#         if sys_Data[i,0]==sys_DataI[index,0]:
#             sys_Data[i,:] = sys_DataI[index,:]
#             index+=1
#     
#     """ Update PQ buses"""                         
#     """Implicit Update (Mismatch): P(T,V)-P_inj,Q(T,V)-Q_inj"""
#     for i in range(n):
#         sys_Data[i,4] = -sys_Data[i,3]
#         sys_Data[i,6] = -sys_Data[i,5]
#         for j in range(n):
#             sys_Data[i,4] += sys_Data[i,1]*sys_Data[j,1]*((sys_G[i,j]*np.cos(sys_Data[i,2]-sys_Data[j,2]))+(sys_B[i,j]*np.sin(sys_Data[i,2]-sys_Data[j,2]))) 
#             sys_Data[i,6] += sys_Data[i,1]*sys_Data[j,1]*((sys_G[i,j]*np.sin(sys_Data[i,2]-sys_Data[j,2]))-(sys_B[i,j]*np.cos(sys_Data[i,2]-sys_Data[j,2]))) 
#             
#     
#     
#     return sys_Data
# 
# def updateE_SysData(sys_Data, sys_G, sys_B, sys_BusType):
#     nodes = sys_BusType.size
#     for i in range(nodes):
#         if sys_BusType[i]=='S':
#             sys_Data[i,3] = 0
#         if sys_BusType[i]!='D':
#             sys_Data[i,5] = 0
#         for j in range(nodes):
#             if sys_BusType[i]=='S':
#                 sys_Data[i,3] += sys_Data[i,1]*sys_Data[j,1]*((sys_G[i,j]*np.cos(sys_Data[i,2]-sys_Data[j,2]))+(sys_B[i,j]*np.sin(sys_Data[i,2]-sys_Data[j,2]))) 
#             if sys_BusType[i]!='D':
#                 sys_Data[i,5] += sys_Data[i,1]*sys_Data[j,1]*((sys_G[i,j]*np.sin(sys_Data[i,2]-sys_Data[j,2]))-(sys_B[i,j]*np.cos(sys_Data[i,2]-sys_Data[j,2])))  
#     return sys_Data
#     
# """
# Code for testing purposes
# """
# """
#  nodes = sys_BusType.size
#     n = nodes-1
#     ng = sys_BusType[sys_BusType!='D'].size
#     nd = sys_BusType[sys_BusType=='D'].size
# 
#     J = np.zeros((n+nd,n+nd))
#     for i in range(n):
#         for j in range(n):                  #(i, j, V_i,             V_j,             T_i,             T_j,             P_i,                               Q_i,                               G_ij,            B_ij)
#             J[i,j] =    Jacobian_PowerFlow_11(i+1, j+1, sys_Data[i+1,1], sys_Data[j+1,1], sys_Data[i+1,2], sys_Data[j+1,2], sys_Data[i+1, 4]+sys_Data[i+1, 3], sys_Data[i+1, 6]+sys_Data[i+1, 5], sys_G[i+1, j+1], sys_B[i+1, j+1])
#     for i in range(n):
#         for j in range(nd): 
#             J[i,j+n] =  Jacobian_PowerFlow_12(i+1, j+ng, sys_Data[i+1,1], sys_Data[j+ng,1], sys_Data[i+1,2], sys_Data[j+ng,2], sys_Data[i+1, 4]+sys_Data[i+1, 3], sys_Data[i+1, 6]+sys_Data[i+1, 5], sys_G[i+1, j+ng], sys_B[i+1, j+ng])
#     for i in range(nd):
#         for j in range(n): 
#             J[i+n,j] =  Jacobian_PowerFlow_21(i+ng, j+1, sys_Data[i+ng,1], sys_Data[j+1,1], sys_Data[i+ng,2], sys_Data[j+1,2], sys_Data[i+ng, 4]+sys_Data[i+ng, 3], sys_Data[i+ng, 6]+sys_Data[i+ng, 5], sys_G[i+ng, j+1], sys_B[i+ng, j+1])
#     for i in range(nd):
#         for j in range(nd): 
#             J[i+n,j+n] =Jacobian_PowerFlow_22(i+ng, j+ng, sys_Data[i+ng,1], sys_Data[j+ng,1], sys_Data[i+ng,2], sys_Data[j+ng,2], sys_Data[i+ng, 4]+sys_Data[i+ng, 3], sys_Data[i+ng, 6]+sys_Data[i+ng, 5], sys_G[i+ng, j+ng], sys_B[i+ng, j+ng])
# """
# """
# test = np.arange(12).reshape(3,4)
# print(test)
# test=np.delete(test, 1,0)
# print(test)
# """
# """
# print(sys_G)
# print()
# print(sys_B)
# """
# 
# """
# dict_BusData = {df_BusData.columns[0]:df_BusData[df_BusData.columns[0]]}
# 
# print(df_BusData['P MW'][df_BusData['Type']=='G'])
# 
# line_From = np.array(df_LineData['From'])
# line_To = np.array(df_LineData['To'])
# busNum = np.array(df_BusData['Bus #'])      
# ind = np.array(df_LineData.columns)
# line_R = np.array(df_LineData[ind[2]])
# """
# 
# """
# ind1 = np.array(line_From==1)
# ind2 = np.array(line_To==2)
# ind3 = ind1 * ind2
# print(-1*line_R[np.array(line_From==1) * np.array(line_To==2)])
# 
# for i in range(df_BusData.shape[0]):
#     print(i)
# #print(df_BusData)
# #print(ind.size)
# #print(build_AdmittanceMatrix(df_LineData))
# """
# """
# test = np.array([[2*i+j for j in range(2)] for i in range(2)], dtype = complex)
# print(test)
# print(test[:,0].size)
# test2 = np.zeros((2,3))
# test2[:,1] = np.array([4,5])
# print(test2)
# #"""
# 
# """
# test1 = 0.5*np.ones((2,2))
# test2 = 1/test1
# #print(test2)
# """
# """
# testshape = np.array([[1,2,3],[1,2,3]])
# print(testshape.shape)
# """
# =============================================================================
