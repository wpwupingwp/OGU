#!/usr/bin/python3
from Bio import SeqIO
from Bio import pairwise2 as p2
import sys

Primer=list()
Out=list()
Unknown=list()
Sum=list()
Sum.append(['all':0])
with open(sys.argv[2],'r') as In:
    Raw=In.read().split(sep='\n')
for line in Raw:
    Primer.append(line.split(sep='\t'))
Primer.pop(-1)
Primer.pop(0)
Sequence=SeqIO.parse(sys.argv[1],'fastq')
for s in Sequence:
    Sum[0][1]+=1
    head=str((s.seq)[0:15])
    Unknown.append(s)
    for p in Primer:
        if head in p[1] or head in p[2]:
                add=[p[0],s]
                Out.append(add)
                Unknown.pop(Unknown.index(s))
                break
for p in Unknown:
    for p in Primer:
        for a in p2.align.localms(head,p[0],1,-1,-0.5,-0.1):
            #1,same -1,different -0.5,gap open -0.1,gap extend
            score=a[2]
        if score>=15:
            add=[p[0],s]
            Out.append(add)
            Unknown.pop(Unknown.index(s))
            break
for cp in Out:
    handle=open(cp[0],'a')
    if cp[0] not in Sum:
        Sum.append([cp[0],1)
    else:
        Sum[Sum.index(cp[0])][1]+=1
    SeqIO.write(cp[1],handle,'fastq')
SeqIO.write(Unknown,'unknown.fastq','fastq')
with open('sum.csv','w') as Out:
    for line in Sum:
        Out.write(line)
