#!/usr/bin/python3

import re
from glob import glob
from Bio import SeqIO

file_list = list(glob('*.fasta'))
pattern = re.compile(r'(\d+);$')
for fasta in file_list:
    new = list()
    for record in SeqIO.parse(fasta, 'fasta'):
        n = re.search(pattern, record.id).group(1)
        n = int(n)
        new.append([n, record])
    # descending
    new.sort(key=lambda x: x[0], reverse=True)
    for index, record in enumerate(new):
        with open('{}.fasta'.format(index+1), 'a') as merge, open(
                '{}.{}'.format(fasta, index+1), 'a') as split:
            SeqIO.write(record[1], merge, 'fasta')
            SeqIO.write(record[1], split, 'fasta')
finish = input('Finish.')
