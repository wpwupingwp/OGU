#!/usr/bin/python3
from Bio.Blast.Applications import NcbiblastnCommandline as nb
from Bio.Blast import NCBIXML as nx 

def runblast():
    cmd=nb(query='unknown.fasta',db='primer',task='blastn-short',evalue=0.001,outfmt=5,out='result')
    stdout,stderr=cmd()
    return 

def parse():
    handle=open('result','r')
    result=list(nx.parse(handle))
    for record in result:
        for item in record.alignments:
            a=item.hsps[0]
            print(a.score,item.hit_def)
            #for hsp in item.hsps:
#                print(item.title,hsp.query,hsp.score)

runblast()
parse()

