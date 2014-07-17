#!/usr/bin/python3
import re
import csv
Value=[]
Line=[]
Row=[]
Fna=input('fna file:\n')
List=input('list file:\n')
with open(Fna,'rb') as In:
    Raw=str(In.read())
with open(List,'r') as List:
    Rows=List.read().split(sep='\n')
if Rows[-1]=='':
    Rows.pop()
Rows.sort()
print(Rows)
Id=re.findall('(?<=\>)[A-Z][a-z]+-cp\d{3}',Raw)
for record in Id:
    Value.append((str(record)).split(sep='-'))
Value.sort()
Line=[]
for n in range(141):
    Line.extend(''.join(['cp','%03d'%(n)]))
Line[0]=''
Out=open('result.csv','wb')
writer=csv.writer(Out)
