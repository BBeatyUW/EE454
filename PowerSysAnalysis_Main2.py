"""
Power Flow Analysis: Support Functions
Split of PQ and PV buses
Created By: 
    Brandon Beaty
    Eliot Nichols
"""
from PowerSysAnalysis_HeaderFile2 import *

"""
########################
    Main Section of Code 
########################    
"""

"""
Loction for Bus and line data
First variable is for file location, Second variable is for sheet name
"""
BusData_Location = ['.\Sample System\system_SampleInput.xlsx', 'BusData'] #Sample Data
LineData_Location = ['.\Sample System\system_SampleInput.xlsx', 'LineData'] #Sample Data
#BusData_Location = ['system_basecase.xlsx', 'BusData'] #BaseCase Data
#LineData_Location = ['system_basecase.xlsx', 'LineData'] #BaseCase Data
tolerance = [.001, .001] #P.U.

"""Data Frame creation for initial Bus and Line Data"""
df_BusData, df_LineData = import_BusAndLineData(BusData_Location, LineData_Location)
nodes = 1+np.arange(df_BusData.shape[0])
sys_G, sys_B = build_AdmittanceMatrix(df_LineData, nodes.size)
sys_LoadP, sys_LoadQ, sys_BusType, sys_PGen, sys_VRef = init_BusData(df_BusData)
sys_Data = init_SysData(sys_LoadP, sys_LoadQ, sys_BusType, sys_PGen, sys_VRef, sys_G, sys_B)
np.set_printoptions(precision=4, suppress=True)
mismatch_P = sys_Data[:,4]
mismatch_Q = sys_Data[:,6]
mismatch_max = [max(abs(mismatch_P)), max(abs(mismatch_Q))]
iteration = 0
while(iteration<100 and mismatch_max>tolerance):
    #print(iteration)
    #print(sys_Data)
    sys_Data = updateI_SysData(sys_Data, sys_G, sys_B, sys_BusType)
    mismatch_P = sys_Data[:,4]
    mismatch_Q = sys_Data[:,6]
    mismatch_max = [max(abs(mismatch_P)), max(abs(mismatch_Q))]
    print(iteration, mismatch_max)
    iteration += 1
sys_Data = updateE_SysData(sys_Data, sys_G, sys_B, sys_BusType)
print(iteration)
print(sys_Data)