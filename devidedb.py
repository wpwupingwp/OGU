#!/usr/bin/python3
from Bio import SeqIO
import sys

area=sys.argv[1][-1]
with open(sys.argv[1],'r') as In:
    List=In.read().split(sep='\n')
List.pop()
handle=open('db'+area,'a')
DB=SeqIO.parse('maindb','fasta')
for record in DB:
    name=record.id.split(sep='_')
    if name[0] in List:
        SeqIO.write(record,handle,'fasta')
