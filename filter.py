#!/usr/bin/python3
from Bio import SeqIO
import sys

raw = list(SeqIO.parse(sys.argv[1], 'fasta'))
filtered = [i for i in raw if len(i.seq) > int(sys.argv[2])]
SeqIO.write(filtered, sys.argv[1]+sys.argv[2], 'fasta')
