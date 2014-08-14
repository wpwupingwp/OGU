#!/usr/bin/python3
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
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
Toblast=list()
for item in Unknown:
    add=SeqRecord(id=item.id,description='',seq=item.seq[:30])
    Toblast.append(add)
SeqIO.write(Toblast,'unknown.fasta','fasta')
with open('sum.csv','w') as Out:
    for key,value in Sum.items():
        Out.write(' '.join([key,str(value),'\n']))
#!/usr/bin/python3
from Bio import SearchIO
from Bio.Blast.Applications import NcbiblastnCommandline as nb
from Bio.Blast import NCBIXML as nx 
import sys

def runblast():
    cmd=nb(query='unknown.fasta',db='primer',task='blastn-short',evalue=0.001,outfmt=5,out='blast.result')
    stdout,stderr=cmd()
    return 

def parse():
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

def newparse():
    global BlastResult
    BlastResult=dict()
    Out=open('blast.log','a')
#    sys.stdout=Out
    results=list(SearchIO.parse('blast.result','blast-xml'))
    for record in results:
        if len(record)==0:
            continue
        else: 
            tophit=record[0]
        print(tophit.id,record.id,'\n',tophit)
        BlastResult[record.id]=tophit.id

runblast()
newparse()

