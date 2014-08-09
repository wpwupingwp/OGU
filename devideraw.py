#!/usr/bin/python3
from Bio import SeqIO
import re

Primer=list()
Out=list()
Unknown=list()
with open('/tmp/work/primer_list.txt','r') as In:
    Raw=In.read().split(sep='\n')
for line in Raw:
    Primer.append(line.split(sep='\t'))
Primer.pop(-1)
Primer.pop(0)
Sequence=SeqIO.parse('/tmp/work/1.fastq','fastq')
n=0
m=0
for s in Sequence:
    n+=1
    head=str((s.seq)[15:30])
    Unknown.append(s)
    for p in Primer:
        if re.search(head,p[1])!=None or re.search(head,p[2])!=None:
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
