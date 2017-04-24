#!/usr/bin/python3

from Bio import SeqIO
from sys import argv

print('Usage:\npython fasta_file list')

fasta_list = list()
with open(argv[2], 'r') as fasta:
    for line in fasta:
        fasta_list.append(line.strip())
handle = open('result.fasta', 'w')
for record in SeqIO.parse(argv[1], 'fasta'):
    info = '>'+record.description 
    if info in fasta_list:
        SeqIO.write(record, handle, 'fasta')
