import numpy as np
import scipy as sc
import math
import openpyxl
import csv
import pandas as pd
import xlrd

BusData_Location = ['.\Sample System\system_SampleInput.xlsx', 'BusData']
LineData_Location = ['.\Sample System\system_SampleInput.xlsx', 'LineData']

def import_BusAndLineData(BusData_Location, LineData_Location):
    BusData = pd.read_excel(BusData_Location[0], sheet_name=BusData_Location[1])
    LineData = pd.read_excel(LineData_Location[0], sheet_name=LineData_Location[1])
    return BusData, LineData


df_BusData, df_LineData = import_BusAndLineData(BusData_Location, LineData_Location)

print(df_BusData)
print(df_LineData)
