#!/usr/bin/python3

import glob
import pandas

#To get all xls filename in present directory
#For example, 'A23-4.xls'
data = list()
name_list = glob.glob('*.xls')
for name in name_list:
#Use first sheet
    sheet = pandas.read_excel(name)
    data.append([name,sheet])
    print(sheet)
