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
#BusData_Sample = ['.\Sample System\system_SampleInput.xlsx', 'BusData'] #Sample Data
#LineData_Sample = ['.\Sample System\system_SampleInput.xlsx', 'LineData'] #Sample Data
#Output_Sample = "SampleData"
BusData_Location = ['system_basecase.xlsx', 'BusData'] #BaseCase Data
LineData_Location_BaseCase = ['system_basecase.xlsx', 'LineData'] #BaseCase Data
LineData_Location_Contingency1 = ['system_basecase.xlsx', 'ContingencyCase1'] #BaseCase Data
LineData_Location_Contingency2 = ['system_basecase.xlsx', 'ContingencyCase2'] #BaseCase Data

tolerance = [.001, .001] #P.U.
S_Base = 100 #MW

PowerFlowAnalysis(BusData_Location, LineData_Location_BaseCase, "Base Case Solution", tolerance, S_Base)
PowerFlowAnalysis(BusData_Location, LineData_Location_Contingency1, "Contingency Case 1 Solution", tolerance, S_Base)
PowerFlowAnalysis(BusData_Location, LineData_Location_Contingency2, "Contingency Case 2 Solution", tolerance, S_Base)