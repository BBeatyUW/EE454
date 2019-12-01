import numpy as np
import scipy as sc
import math
import openpyxl
import csv
import pandas as pd
import xlrd

df = pd.read_excel(r'system_SampleInput.xlsx',sheet_name='BusData')

print(df)
