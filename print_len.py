#!/usr/bin/python3

from Bio import SeqIO
from glob import glob


for fasta in glob('*.fa*'):
    print(fasta)
    for record in SeqIO.parse(fasta, 'fasta'):
        print(record.id, len(record))
    print()
