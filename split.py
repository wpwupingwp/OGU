#!/usr/bin/python

from glob import glob
from Bio import SeqIO

filelist = list(glob('*.fasta'))
for fasta in filelist:
    for record in SeqIO.parse(fasta, 'fasta'):
        if record.id.count('BOP') == 1:
            with open(fasta+'.bop', 'a') as output:
                SeqIO.write(record, output, 'fasta')
        else:
            with open(fasta+'.nobop', 'a') as output:
                SeqIO.write(record, output, 'fasta')
