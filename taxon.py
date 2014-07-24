#!/usr/bin/python3

import csv
Id=dict()
Rank=dict()
Name=dict()
List=[]
Data=[]
with open('./test/id','r') as In:
    Raw=list(In.readlines())
    for record in Raw:
        add=record.replace('\n','').split(sep=' ')
        Id[add[0]]=add[1]
with open('./test/rank','r') as In:
    Raw=list(In.readlines())
    for record in Raw:
        add=record.replace('\n','').split(sep=' ')
        Rank[add[0]]=add[1]
        if add[1]=='species':
            List.append(add[0])
with open('./test/name','r') as In:
    Raw=list(In.readlines())
    for record in Raw:
        add=record.replace('\n','').split(sep='|')
        if add[0] not in Name or add[3]=='scientific name':
            Name[add[0]]=add[1]
for species in List:
    record=[species,]
    while Id[species]!='1' :
#        if Rank[Id[species]] in ['species','genus','family','order','class','phylum','kingdom'] :
        record.append(Id[species])
        species=Id[species]
    if '33090' in record:
        record.pop()
        record.pop()
        Data.append(record[::-1])

writer=csv.writer(open('out.csv','w',newline=''))
writer2=csv.writer(open('out2.csv','w',newline=''))
for item in Data:
    writer.writerow(item)
    item2=[]
    for word in item:
        word=Name[word]
        item2.append(word)
    writer2.writerow(item2)
        
