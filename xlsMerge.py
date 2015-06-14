#!/usr/bin/python3

import glob
import xlrd

#To get all xls filename in present directory
#For example, 'A23-4.xls'
data = list()
name_list = glob.glob('*.xls')
for name in name_list:
    book = xlrd.open_workbook(name)
#Use first sheet
    sheet = book.sheet_by_index(0)
    data.append([name,sheet])
print(data)
