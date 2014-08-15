#!/usr/bin/python3
import sys
import re
from Bio import SeqIO

with open(sys.argv[1],'r') as fna:
    In=fna.read()
with open(sys.argv[2],'r') as qual:
    Out=qual.read()
Rawinfo=re.findall('(?<=\>)[0-9a-zA-Z_].*[0-9]{5}',In)
for name in Rawinfo:
    info=name[-15:]
    Out=re.sub(info,name,Out,count=1)
with open(sys.argv[2],'w') as out:
    out.write(Out)
fna=list(SeqIO.parse(sys.argv[1],'fasta'))
qual=list(SeqIO.parse(sys.argv[2],'qual'))
fna.sort(key=lambda x:x.id)
qual.sort(key=lambda x:x.id)
SeqIO.write(fna,sys.argv[1],'fasta')
SeqIO.write(qual,sys.argv[2],'qual')
