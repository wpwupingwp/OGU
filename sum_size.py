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
        print(seqs)
        new.append([seqs, record])
    # descending
    new.sort(key=lambda x: x[0], reverse=True)

    for m, n in enumerate(new):
        n[0] = m + 1
    for i in new:
        print(i)
    for n in range(1, len(new)):
        with open('{}.fasta'.format(n), 'a') as merge, open(
                '{}.{}'.format(fasta, n), 'a') as split:
            SeqIO.write(new[n][1], merge, 'fasta')
            SeqIO.write(new[n][1], split, 'fasta')
finish = input('Finish.')
