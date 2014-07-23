#!/usr/bin/python3

Id=dict()
Rank=dict()
List=[]
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
            List.append(add)
print(List)
