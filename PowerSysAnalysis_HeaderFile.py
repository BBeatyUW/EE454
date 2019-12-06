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
                sys_Y[i][j] = -np.sum(line_Y[np.multiply(np.array(line_From==i+1), np.array(line_To==j+1))])
            else: #i>j =[j][i]
               sys_Y[i][j] = sys_Y[j][i] 
    sys_G = sys_Y.real
    sys_B = sys_Y.imag
    return sys_Y, sys_G, sys_B

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



"""
Returns Appearent, real, and reactive power flowing a line.
"""
def PowerFlow (V_i,T_i,V_j,T_j,B_tot,y_ij):
    I_ij = y_ij * (V_i * np.cos(T_i) + 1j * V_i * np.sin(T_i) - V_j * np.cos(T_j) - 1j * V_j
                   * np.sin(T_j)) + (1j*B_tot / 2) * (V_i * np.cos(T_i) + 1j * V_i * np.sin(T_i))
    S_ij = (V_i*np.cos(T_i)+1j*V_i*np.sin(T_i)) * (I_ij.conjugate())
    return abs(S_ij), S_ij.real, S_ij.imag


"""
This fuction collects all the needed line data and turns it into lists for exporting.
LD_val is a matrix made from the dataframe that is given.
S_ij is a list that stores all the apparent power flowing in the line (s_ij)
P_ij is a list that stores all the apparent power flowing in the line (p_ij)
Q_ij is a list that stores all the apparent power flowing in the line (q_ij)
i_buses is a list of the buses that are used as a starting point
j_buses is a list of the buses that are used as an ending point
violataion list stores whether the line power limit was violated or not
"""
def LineFlowResults (sys_Data, LineData, sys_Y):
    LD_val = LineData.values
    S_ij = []
    P_ij = []
    Q_ij = []
    i_buses = LD_val[0:,0]
    j_buses = LD_val[0:,1]
    violation = []
    for i in range (0, len(LineData)):
        B_tot = (LD_val[i:i+1,4])
        V_i = sys_Data[(int(LD_val[i:i+1,0]))-1:(int(LD_val[i:i+1,0])),1]
        T_i = sys_Data[(int(LD_val[i:i + 1, 0])) - 1:(int(LD_val[i:i + 1, 0])), 2]
        V_j = sys_Data[(int(LD_val[i:i + 1, 1])) - 1:(int(LD_val[i:i + 1, 1])), 1]
        T_j = sys_Data[(int(LD_val[i:i + 1, 1])) - 1:(int(LD_val[i:i + 1, 1])), 2]
        y_ij = -1 * ( sys_Y[(int(LD_val[i:i + 1, 0]) - 1):(int(LD_val[i:i + 1, 0])),
               (int(LD_val[i:i + 1, 1]) - 1)])

        PowerFlow(V_i, T_i, V_j, T_j, B_tot, y_ij)
        s_ij, p_ij, q_ij = PowerFlow(V_i,T_i,V_j,T_j,B_tot,y_ij)
        
        print('Line Violation: ', s_ij, LD_val[i:i+1,5])
        if s_ij*100 < (LD_val[i:i+1,5]):
            violation.append('FALSE')
        else:
            violation.append('TRUE')
        S_ij.append(100 * float(s_ij))
        P_ij.append(100 * float(p_ij))
        Q_ij.append(100 * float(q_ij))
    return S_ij, P_ij, Q_ij, violation, i_buses, j_buses

"""
This fuction collects all the needed bus data and turns it into lists for exporting.
"""
def BusResults(sys_Data):
    V_violate = []
    for i in range (0,len(sys_Data)):
        if  (sys_Data[i:i+1,1] <= 1.05 and sys_Data[i:i+1,1] >= 0.95):
            V_violate.append('FALSE')
        else:
            V_violate.append('TRUE')
    bus_nums = sys_Data[0:,0]
    bus_v = sys_Data[0:,1]
    bus_deg = 180* sys_Data[0:,2] / np.pi
    bus_p = 100 * sys_Data[0:,3]
    bus_q = 100 * sys_Data[0:,5]
    return bus_nums.astype(int), bus_v, bus_deg, bus_p, bus_q, V_violate

"""
This function collects all the data obtained and outputs it to a single excel file with multiple sheets.
"""
def DataOutput(sys_Data, LineData, sys_Y,iteration_list,mismatch_P_list,mismatch_Q_list,max_P_bus,max_Q_bus):

    y_matrix = pd.DataFrame(data=sys_Y)
    g_matrix = pd.DataFrame(data=sys_Y.real)
    b_matrix = pd.DataFrame(data=sys_Y.imag)

    bus_nums, bus_v, bus_deg, bus_p, bus_q, V_violate = BusResults(sys_Data)
    S_ij, P_ij, Q_ij, S_violation, i_buses, j_buses = LineFlowResults(sys_Data, LineData, sys_Y)

    df_Convergence = pd.DataFrame({'Iteration': iteration_list, 'Max P Misr': mismatch_P_list, 'P Bus': max_P_bus
                                    , 'Max Q Misr': mismatch_Q_list, 'Q Bus': max_Q_bus})

    df_BusOutput = pd.DataFrame({'Bus Number':bus_nums,'V (pu)':bus_v,'Angle (deg)':bus_deg
                                    ,'P inj. (MW)':bus_p,'Q inj. (MVar)':bus_q,
                                 'Voltage Violation':V_violate})

    df_LineOutput = pd.DataFrame({'From Bus': i_buses, 'To Bus': j_buses, 'P Flow (MW)': P_ij, 'Q Flow (MVar)': Q_ij, 'S Flow (MVA)': S_ij,
                                 'Line MVA Violation': S_violation})

    writer = pd.ExcelWriter("PowerSystemData.xlsx", engine='xlsxwriter')

    df_BusOutput.to_excel(writer, sheet_name='BusData',startrow=1,header=False,index=False)
    df_LineOutput.to_excel(writer, sheet_name='LineData', startrow=1, header=False, index=False)
    df_Convergence.to_excel(writer, sheet_name='ConvergenceHistory', startrow=1, header=False, index=False)
    y_matrix.to_excel(writer, sheet_name='Y_matrix(Admittance)', startrow=0,header=False,index=False)
    g_matrix.to_excel(writer, sheet_name='G_matrix', startrow=0,header=False,index=False)
    b_matrix.to_excel(writer, sheet_name='B_matrix', startrow=0,header=False,index=False)

    workbook = writer.book
    busworksheet = writer.sheets['BusData']
    lineworksheet = writer.sheets['LineData']
    convergencesheet = writer.sheets['ConvergenceHistory']

    header_format = workbook.add_format({'bold': True,
    'text_wrap': True,
    'valign': 'top',
    'fg_color': '#D7E4BC',
    'border': 1})

    for col_num, value in enumerate(df_BusOutput.columns.values):
        busworksheet.write(0, col_num,value,header_format)
    for col_num, value in enumerate(df_LineOutput.columns.values):
        lineworksheet.write(0, col_num, value, header_format)
    for col_num, value in enumerate(df_Convergence.columns.values):
        convergencesheet.write(0, col_num, value, header_format)

    writer.save()
