#!/usr/bin/python3
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline as nb
import sys

def RunBlast():
    cmd=nb(query='unknown.fasta',db='primer',task='blastn-short',evalue=0.001,outfmt=5,out='blast.result')
    stdout,stderr=cmd()
    return 

def OldParse():
    global BlastResult
    BlastResult=dict()
    In=open('blast.result','r')
    Out=open('blast.log','a')
    sys.stdout=Out
    results=list(nx.parse(In))
    for record in results:
        for hit in record.alignments:
            a=hit.hsps[0]
            if a.score<15:
                continue
            BlastResult[record.query]=hit.hit_def
            print(record.query,hit.hit_def,'\n',a)

def NewParse():
    Out=open('blast.log','a')
    sys.stdout=Out
    results=list(SearchIO.parse('blast.result','blast-xml'))
    for record in results:
        if len(record)==0:
            continue
        else: 
            tophit=record[0]
            #ignore multiple hsps
        print(tophit.id,record.id,'\n',tophit)
        BlastResult[record.id]=tophit.id

#main
Out=list()
Unknown=list()
Toblast=list()
BlastResult=dict()
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
for item in Unknown:
    add=SeqRecord(id=item.id,description='',seq=item.seq[:30])
    Toblast.append(add)
SeqIO.write(Unknown,'unknown.fasta','fasta')
#SeqIO.write(Toblast,'unknown.fasta','fasta')
RunBlast()
NewParse()
for index,record in enumerate(Unknown):
    if record.id in BlastResult:
        Unknown.pop(index)
for cp in Out:
    handle=open(''.join([cp[0],'.fastq']),'a')
    Sum[cp[0]]+=1
    SeqIO.write(cp[1],handle,'fastq')
Sum['unknown']=len(Unknown)
Sum['blasted']=len(BlastResult)
Sum['all']=all
with open('sum.csv','w') as Out:
    for key,value in Sum.items():
        Out.write(' '.join([key,str(value),'\n']))
