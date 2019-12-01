import numpy as np
#import scipy as sc
#import math
#import openpyxl
#import csv
import pandas as pd
#import xlrd

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

"""
def init_BusData(BusData):
    col = np.array(BusData.columns)
    sys_LoadP = np.array(BusData[col[1]])
    sys_LoadQ = np.array(BusData[col[2]])
    sys_Type = np.array(BusData[col[3]])
    sys_P = np.array(BusData[col[4]])
    sys_V = np.array(BusData[col[5]])
    return sys_BusData

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

"""Data Frame creation for initial Bus and Line Data"""
df_BusData, df_LineData = import_BusAndLineData(BusData_Location, LineData_Location)
sys_G, sys_B = build_AdmittanceMatrix(df_LineData, df_BusData.shape[0])

"""
print(sys_G)
print()
print(sys_B)
"""

"""
Code for testing purposes
"""
line_From = np.array(df_LineData['From'])
line_To = np.array(df_LineData['To'])
busNum = np.array(df_BusData['Bus #'])      
ind = np.array(df_LineData.columns)
line_R = np.array(df_LineData[ind[2]])

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

test = np.array([[i for j in range(2)] for i in range(2)], dtype = complex)
print(test)
print(test[:,0])

test1 = 0.5*np.ones((2,2))
test2 = 1/test1
#print(test2)
"""