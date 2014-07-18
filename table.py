#!/usr/bin/python3
import re
import csv
Value=[]
Line=[]
Row=[]
Writeline=[]
Out=[]
Fna=input('fna file:\n')
List=input('list file:\n')
with open(Fna,'rb') as In:
    Raw=str(In.read())
with open(List,'r') as List:
    Rows=List.read().split(sep='\n')
Rows.sort()
Rows.pop(0)
for n in range(141):
    Line.append(''.join(['cp','%03d'%(n)]))
Line[0]=''
Out.append(Line)
fill=['0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0']
for n in range(len(Rows)):
    if Rows[n]!='':
        Add=[Rows[n],]
        Add.extend(fill)
        Out.append(Add)
Id=re.findall('(?<=\>)[A-Z][a-z]+-cp\d{3}',Raw)
for record in Id:
    Value.append((str(record)).split(sep='-'))
Value.sort()
for item in Value:
    if item[0] in Rows:
        y=Rows.index(item[0])+0
        x=Line.index(item[1])
        Out[y][x]=1
handle=open('output.csv','w')
writer=csv.writer(handle)
for line in Out:
    writer.writerow(line)

