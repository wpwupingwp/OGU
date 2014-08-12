#!/usr/bin/python3
from Bio import SeqIO
from Bio import pairwise2 as p2
import sys

Primer=list()
Out=list()
Unknown=list()
with open(sys.argv[2],'r') as In:
    Raw=In.read().split(sep='\n')
for line in Raw:
    Primer.append(line.split(sep='\t'))
Primer.pop(-1)
Primer.pop(0)
Sequence=SeqIO.parse(sys.argv[1],'fastq')
n=0
m=0
for s in Sequence:
    n+=1
    head=str((s.seq)[0:20])
    Unknown.append(s)
    for p in Primer:
        for a in p2.align.localmx(head,p[0],1,-1,-0.5,-0.1):
            score=a[2]
            if score>=15:
                add=[p[0],s]
                Out.append(add)
                Unknown.pop(Unknown.index(s))
                m+=1
                break
print(n,m,len(Unknown))
#for cp in Out:
#    handle=open(cp[0],'a')
#    SeqIO.write(cp[1],handle,'fastq')
#SeqIO.write(Unknown,'unknown.fastq','fastq')
