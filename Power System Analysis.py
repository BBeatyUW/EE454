import numpy as np
from numpy.linalg import inv
import pandas as pd

"""
Supporting Fucntions
TODO:
    Section these off into a seperate file
    and import them into main file
"""

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
[3] - Reactive Power (Q)
"""
def init_SysData(sys_LoadP, sys_LoadQ, sys_BusType, sys_PGen, sys_VRef):
    n= sys_LoadP.size
    sys_Data = np.zeros((n,4))
    sys_Data[:,0] = sys_VRef.copy()
    sys_Data[:,1] = np.zeros(n)
    """
    for i in range(n):
        for j in range(n):
            if sys_BusType[i]=='G':
                sys_Data[i,2]=sys_PGen[i]
            else:
                sys_Data[i, 2] += sys_Data[i, 0]*sys_Data[j,0]*(sys_G[i,j]*np.cos(sys_Data[i,1]-sys_Data[j,1])+sys_B[i,j]*np.sin(sys_Data[i,1]-sys_Data[j,1]))
            sys_Data[i, 3] += sys_Data[i, 0]*sys_Data[j,0]*(sys_G[i,j]*np.sin(sys_Data[i,1]-sys_Data[j,1])-sys_B[i,j]*np.cos(sys_Data[i,1]-sys_Data[j,1]))
    """
    sys_Data[:,2] = -1*sys_LoadP[sys_BusType=='D']
    sys_Data[:,3] = -1*sys_LoadQ[sys_BusType=='D'] 
       
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
            J = V_i*V_j*(G_ij*np.sin(T_i-T_j)-B_ij*np.cos(T_i-T_j))
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
[2] - Active Power (P)
[3] - Reactive Power (Q),
As well as, the systems G and B matrices, and node types
Returns the updated array
"""
def update_SysData(sys_Data, sys_G, sys_B, sys_BusType, sys_LoadP, sys_LoadQ):
    n=sys_Data.shape[0]
    """ Determine Jacobian """
    J = np.zeros((2*(n-1),2*(n-1)))
    for i in range(2*(n-1)):
        for j in range(2*(n-1)):
            J[i,j] = Jacobian_PowerFlow(i, j, (n-1), sys_Data[(i+1)%n,0], sys_Data[(j+1)%n,0],\
                 sys_Data[(i+1)%n,1], sys_Data[(j+1)%n,1], sys_Data[(i+1)%n,2], sys_Data[(i+1)%n,3],\
                 sys_G[(i+1)%n,(j+1)%n], sys_B[(i+1)%n,(j+1)%n], )
    
    """ Determine inverse of Jacobian """
    J_inv = inv(J)
    
    """ Calculate Delta V's, Theta's, and update """
    
    Delta = J_inv @ np.append(sys_Data[1:n,2], sys_Data[1:n,3])       
    for i in range((n-1)):
        sys_Data[i+1, 0] += Delta[i+(n-1)]
        sys_Data[i+1, 1] += Delta[i]
        
    """ Update P,Q """
    for i in range(n):
        sys_Data[i, 2] = -sys_LoadP[i]
        sys_Data[i, 3] = -sys_LoadQ[i]
        for j in range(n):
            sys_Data[i, 2] += sys_Data[i, 0]*sys_Data[j,0]*(sys_G[i,j]*np.cos(sys_Data[i,1]-sys_Data[j,1])+sys_B[i,j]*np.sin(sys_Data[i,1]-sys_Data[j,1]))
            sys_Data[i, 3] += sys_Data[i, 0]*sys_Data[j,0]*(sys_G[i,j]*np.sin(sys_Data[i,1]-sys_Data[j,1])-sys_B[i,j]*np.cos(sys_Data[i,1]-sys_Data[j,1])) 
        
    
    return sys_Data

"""
########################
    Main Section of Code 
########################    
"""

"""
Loction for Bus and line data
First variable is for file location, Second variable is for sheet name
"""
BusData_Location = ['.\Sample System\system_SampleInput.xlsx', 'BusData']
LineData_Location = ['.\Sample System\system_SampleInput.xlsx', 'LineData']
tolerance = [.1, .1]

"""Data Frame creation for initial Bus and Line Data"""
df_BusData, df_LineData = import_BusAndLineData(BusData_Location, LineData_Location)
nodes = 1+np.arange(df_BusData.shape[0])
sys_G, sys_B = build_AdmittanceMatrix(df_LineData, nodes.size)
sys_LoadP, sys_LoadQ, sys_BusType, sys_PGen, sys_VRef = init_BusData(df_BusData)
sys_Data = init_SysData(sys_LoadP, sys_LoadQ, sys_BusType, sys_PGen, sys_VRef)

mismatch_P = np.max(sys_LoadP-sys_Data[:,2])
mismatch_Q = np.max(sys_LoadQ-sys_Data[:,3])
mismatch = [abs(mismatch_P), abs(mismatch_Q)]
iteration = 0
convergence = [['Iteration', 'Max P Mismatch', 'P Bus', 'Max Q Mismatch', 'Q Bus'],\
               [iteration, mismatch[0], nodes[sys_LoadP==mismatch_P+sys_Data[:,2]],\
               mismatch[1], nodes[sys_LoadQ==mismatch_Q+sys_Data[:,3]]]]

while(mismatch>tolerance and iteration < 5):
    print([iteration, mismatch[0], nodes[sys_LoadP==mismatch_P+sys_Data[:,2]],\
               mismatch[1], nodes[sys_LoadQ==mismatch_Q+sys_Data[:,3]]])
    update_SysData(sys_Data, sys_G, sys_B, sys_BusType, sys_LoadP, sys_LoadQ)
    mismatch_P = np.max(sys_LoadP-sys_Data[:,2])
    mismatch_Q = np.max(sys_LoadQ-sys_Data[:,3])
    mismatch = [abs(mismatch_P), abs(mismatch_Q)]
    iteration+=1
    convergence.append([iteration, mismatch[0], nodes[sys_LoadP==mismatch_P+sys_Data[:,2]],\
               mismatch[1], nodes[sys_LoadQ==mismatch_Q+sys_Data[:,3]]])



"""
print(sys_G)
print()
print(sys_B)
"""
"""
Code for testing purposes
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