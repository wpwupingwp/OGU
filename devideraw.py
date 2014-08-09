#!/usr/bin/python3
from Bio import SeqIO

Primer=list()
with open('/tmp/work/primer_list.txt','r') as In:
    Raw=In.read().split(sep='\n')
for line in Raw:
    Primer.append(line.split(sep='\t'))
Primer.pop(-1)
Primer.pop(0)
Sequence=SeqIO.parse('/tmp/work/1.fastq','fastq')
for s in Sequence:
    head=str((s.seq)[15:35])
    for p in Primer:
        if head in p[1] or head in p[2]:
            print(head,p,'\n')
