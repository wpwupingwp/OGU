#!/usr/bin/python3

from Bio import SeqIO
import re

Primer=list()
with open('primer_list.txt','r') as In:
    Raw=In.readlines()
for line in Raw:
    Primer.append(line.split(sep='\t'))
print(Primer)
#Sequence=SeqIO.parse('1.fastq','fastq')


