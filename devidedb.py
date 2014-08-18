#!/usr/bin/python3
from Bio import SeqIO
import sys
import re

area=re.search('\d+',sys.argv[2]).group()
with open(sys.argv[2],'r') as In:
    List=In.read().split(sep='\n')
List.pop()
handle=open('db'+area,'a')
DB=SeqIO.parse(sys.argv[1],'fasta')
for record in DB:
    name=record.id.split(sep='_')
    if name[0] in List:
        SeqIO.write(record,handle,'fasta')
