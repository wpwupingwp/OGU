#!/usr/bin/python3

import csv
Id=dict()
Rank=dict()
List=[]
Data=[]
with open('id','r') as In:
    Raw=list(In.readlines())
    for record in Raw:
        add=record.replace('\n','').split(sep=' ')
        Id[add[0]]=add[1]
with open('rank','r') as In:
    Raw=list(In.readlines())
    for record in Raw:
        add=record.replace('\n','').split(sep=' ')
        Rank[add[0]]=add[1]
        if add[1]=='species':
            List.append(add[0])
for species in List:
    record=[species,]
    while Id[species]!='1':
        record.append(Id[species])
        species=Id[species]
    Data.append(record)

writer=csv.writer(open('out.csv','w',newline=''))
for item in Data:
    if '33090' in item:
        writer.writerow(item)
