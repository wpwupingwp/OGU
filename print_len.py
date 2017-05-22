#!/usr/bin/python3

from Bio import SeqIO
from glob import glob


for fasta in glob('*.fa*'):
    print(fasta)
    for record in SeqIO.parse(fasta, 'fasta'):
        length = len(record)
        if length > 10000:
            print(record.id, '\t', length)
    print()
