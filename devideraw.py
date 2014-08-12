#!/usr/bin/python3
from Bio import SeqIO
from Bio import pairwise2 as p2
import sys
import re

Primer=list()
Out=list()
Unknown=list()
Sum={'cp{:03d}'.format(n+1):0 for n in range(140)}
with open(sys.argv[2],'r') as In:
    Raw=In.read().split(sep='\n')
for line in Raw:
    Primer.append(line.split(sep='\t'))
Primer.pop(-1)
Primer.pop(0)
Unknown=list(SeqIO.parse(sys.argv[1],'fastq'))
all=len(Unknown)
for index,record in enumerate(Unknown):
    head=str((record.seq)[0:15])
    for p in Primer:
        score=0
        if re.search(head,p[1])!=None or re.search(head,p[2])!=None:
            add=[p[0],record]
            Out.append(add)
            Unknown.pop(index)
            break
        else:
            for a in p2.align.localms(head,p[0],1,-1,-0.5,-0.1):
                #1,same -1,different -0.5,gap open -0.1,gap extend
                score=a[2]
            if score>=15:
                add=[p[0],record]
                Out.append(add)
                Unknown.pop(index)
                break
for cp in Out:
    handle=open(''.join([cp[0],'.fastq']),'a')
    Sum[cp[0]]+=1
    print(Sum[cp[0]])
    SeqIO.write(cp[1],handle,'fastq')
Sum['unknown']=len(Unknown)
Sum['all']=all
SeqIO.write(Unknown,'unknown.fastq','fastq')
with open('sum.csv','w') as Out:
    for key,value in Sum.items():
        Out.write(' '.join([key,str(value),'\n']))
