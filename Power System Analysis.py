import numpy as np
import scipy as sc
import math
import openpyxl
import csv
import pandas as pd
import xlrd

# Loction for Bus and line data
# First variable is for file location, Second variable is for sheet name
BusData_Location = ['.\Sample System\system_SampleInput.xlsx', 'BusData']
LineData_Location = ['.\Sample System\system_SampleInput.xlsx', 'LineData']

# Imports Bus and line data from excel sheets
# Takes in an array containing ['File Location', 'Sheet Name']
# Returns two panda data frames for the bus and line data
def import_BusAndLineData(BusData_Location, LineData_Location):
    BusData = pd.read_excel(BusData_Location[0], sheet_name=BusData_Location[1])
    LineData = pd.read_excel(LineData_Location[0], sheet_name=LineData_Location[1])
    return BusData, LineData

# Parses LineData data frame for use in a building admittance matrix
# Test function of data parsing
# Returns line_:
#   Num - number of lines
#   From - line origin
#   To - line destination
#   R - line resistance
#   X - line reactance
#   B - shunt capactiance
#   Fmax - maximum power for the line
def parse_LineData(LineData):
    col = np.array(LineData.columns)
    line_From = np.array(LineData[col[0]])
    line_To = np.array(LineData[col[1]])
    line_R = np.array(LineData[col[2]])
    line_X = np.array(LineData[col[3]])
    line_B = np.array(LineData[col[4]])
    line_Fmax = np.array(LineData[col[5]])
    line_Num = line_From.size
    return line_Num, line_From, line_To, line_R, line_X, line_B, line_Fmax

#Builds G and B matrices to be used in Power Flow calculations
def build_AdmittanceMatrix(LineData):
    col = np.array(LineData.columns)
    line_From = np.array(LineData[col[0]])
    line_To = np.array(LineData[col[1]])
    line_R = np.array(LineData[col[2]])
    line_X = np.array(LineData[col[3]])
    line_B = np.array(LineData[col[4]])
    line_Fmax = np.array(LineData[col[5]])
    line_Num = line_From.size
    
    
    return sys_G, sys_B

#Main Section of code
df_BusData, df_LineData = import_BusAndLineData(BusData_Location, LineData_Location)



#Code for testing purposes
line_Bus = np.array(df_LineData['From'])
busNum = np.array(df_BusData['Bus #'])      
ind = np.array(df_LineData.columns)

#print(df_BusData)
#print(ind.size)
print(parse_LineData(df_LineData))