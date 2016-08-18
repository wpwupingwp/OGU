#!/usr/bin/python3

import sys
from Bio import SeqIO

wanted = ['Ta1a', 'Ta2a', 'Ta3a', 'Ta4a', 'Ta5a', 'Ta6a', 'Ta7a']
present = SeqIO.parse(sys.argv[1], 'fasta')
handle = open('filtered.fasta', 'a')
for sequence in present:
    for i in wanted:
        if i in sequence.id:
            SeqIO.write(sequence, handle, 'fasta')
handle.close()
