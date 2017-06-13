#!/usr/bin/python3

import re
from glob import glob
from Bio import SeqIO

file_list = list(glob('*.fasta'))
pattern = re.compile(r'(\d+);$')
for fasta in file_list:
    new = list()
    for record in SeqIO.parse(fasta, 'fasta'):
        seqs = re.search(pattern, record.id).group(1)
        seqs = int(seqs)
        new.append([seqs, record])
    # descending
    new.sort(key=lambda x: x[0], reverse=True)

    for m, n in enumerate(new):
        n[0] = m + 1
    for n in range(len(new)):
        with open('{}.fasta'.format(n+1), 'a') as merge, open(
                '{}.{}'.format(fasta, n+1), 'a') as split:
            SeqIO.write(new[n][1], merge, 'fasta')
            SeqIO.write(new[n][1], split, 'fasta')
finish = input('Finish.')
