"""
Power Flow Analysis: Main Function
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
#BusData_Location = ['.\Sample System\system_SampleInput.xlsx', 'BusData'] #Sample Data
#LineData_Location = ['.\Sample System\system_SampleInput.xlsx', 'LineData'] #Sample Data
BusData_Location = ['system_SampleInput.xlsx', 'BusData'] #BaseCase Data
LineData_Location = ['system_SampleInput.xlsx', 'LineData'] #BaseCase Data
tolerance = [.001, .001] #P.U.
S_Base = 100 #MW

"""Data Frame creation for initial Bus and Line Data"""
df_BusData, df_LineData = import_BusAndLineData(BusData_Location, LineData_Location)
nodes = 1+np.arange(df_BusData.shape[0])
sys_Y, sys_G, sys_B = build_AdmittanceMatrix(df_LineData, nodes.size)
sys_BusNum, sys_LoadP, sys_LoadQ, sys_BusType, sys_PGen, sys_VRef = init_BusData(df_BusData)
sys_Data = init_SysData(sys_BusNum, sys_LoadP, sys_LoadQ, sys_BusType, sys_PGen, sys_VRef, sys_G, sys_B, S_Base)
np.set_printoptions(precision=4, suppress=True)
mismatch_P = sys_Data[:,4]
mismatch_Q = sys_Data[:,6]
mismatch_max = [max(abs(mismatch_P)), max(abs(mismatch_Q))]
iteration = 0
iteration_list = []
mismatch_P_list = []
mismatch_Q_list = []
max_P_bus = []
max_Q_bus = []
while(iteration<15 and mismatch_max>tolerance):
    print(iteration)
    print(sys_Data)


    """
    This chunk below is needed for output. Don't delete, only move
    """
    iteration_list.append(iteration)
    sys_Data = update_SysData(sys_Data, sys_G, sys_B, sys_BusType)
    mismatch_P = sys_Data[:,4]
    mismatch_Q = sys_Data[:,6]
    bus_P, = np.where(mismatch_P == max(abs(mismatch_P)))
    if len(bus_P) == 0:
        bus_P, = np.where(mismatch_P == -1*max(abs(mismatch_P)))
    max_P_bus.append(int(bus_P+1))
    bus_Q, = np.where(mismatch_Q == max(abs(mismatch_Q)))
    if len(bus_Q) == 0:
        bus_Q, = np.where(mismatch_Q == -1*max(abs(mismatch_Q)))
    max_Q_bus.append(int(bus_Q+1))
    mismatch_P_list.append(max(abs(mismatch_P)))
    mismatch_Q_list.append(max(abs(mismatch_Q)))


    mismatch_max = [max(abs(mismatch_P)), max(abs(mismatch_Q))]
    print(iteration, mismatch_max)


    """This line is also needed for output. Don't delete, only move"""
    iteration += 1



print(iteration)
print(sys_Data)

"""This line exports final data to excel. Don't delete, but feel free to move."""
DataOutput(sys_Data, df_LineData, sys_Y,iteration_list,mismatch_P_list,mismatch_Q_list,max_P_bus,max_Q_bus)
