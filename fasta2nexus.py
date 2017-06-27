#!/usr/bin/python3

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from sys import argv

SeqIO.convert(argv[1], 'fasta', argv[2], 'nexus', alphabet=IUPAC.ambiguous_dna)
