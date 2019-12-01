"""
Power Flow Analysis: Support Functions
Created By: 
    Brandon Beaty
    Eliot Nichols
"""
from PowerSysAnalysis_HeaderFile import *

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
               [iteration, mismatch[0], (nodes[sys_LoadP==mismatch_P+sys_Data[:,2]])[0],\
               mismatch[1], (nodes[sys_LoadQ==mismatch_Q+sys_Data[:,3]])[0]]]

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