#!/usr/bin/python3

import sys
from Bio import SeqIO

length = list()
contigs = SeqIO.parse(sys.argv[1], 'fasta')
for contig in contigs:
    length.append(len(contig.seq))
length.sort()

total = 0
for i in length:
    print(i)
    total += i
print('total:\t',  total)
