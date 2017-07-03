#!/usr/bin/python3

from Bio import SeqIO
from sys import argv

raw = list(SeqIO.parse(argv[1], 'fasta'))
longest = max(raw, key=lambda x: len(x.seq))
SeqIO.write(longest, argv[1]+'.longest', 'fasta')
