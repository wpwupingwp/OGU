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
    sys.stdout=Out
    results=list(SearchIO.parse('blast.result','blast-xml'))
    for record in results:
        try:
            print(record[0])
        except:
            pass

runblast()
newparse()

