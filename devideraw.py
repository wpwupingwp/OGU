#!/usr/bin/python3
from Bio import SeqIO
from Bio import pairwise2 as p2
from Bio.pairwise2 import format_alignment as fa
import sys

def pairwise():
    #1,same -1,different -0.5,gap open -0.1,gap extend
    for index,record in enumerate(Unknown):
        if len(record.seq)<100:
            continue
        head=str((record.seq)[2:17])
        for p in Primer:
            aln=p2.align.localms(head,str(p.seq),1,-1,-0.5,-0.1)   #bug 
            score=aln[0][2]
            if score>=15:
                add=[p.id,record]
                Out.append(add)
                Unknown.pop(index)
                break

#main
Out=list()
Unknown=list()
Sum={'cp{:03d}'.format(n+1):0 for n in range(140)}
Primer=list(SeqIO.parse(sys.argv[2],'fasta'))
Unknown=list(SeqIO.parse(sys.argv[1],'fastq'))
all=len(Unknown)
for index,record in enumerate(Unknown):
    head=str((record.seq)[2:17])
    for p in Primer:
        score=0
        if head in p.seq:
            add=[p.id[:-1],record]
            Out.append(add)
            Unknown.pop(index)
            break
#        else:
#            pairwise()
for cp in Out:
    handle=open(''.join([cp[0],'.fastq']),'a')
    Sum[cp[0]]+=1
    SeqIO.write(cp[1],handle,'fastq')
Sum['unknown']=len(Unknown)
Sum['all']=all
SeqIO.write(Unknown,'unknown.fastq','fastq')
with open('sum.csv','w') as Out:
    for key,value in Sum.items():
        Out.write(' '.join([key,str(value),'\n']))
