#!/usr/bin/python3

from Bio import SeqIO
from Bio import Alphabet
import sys

SeqIO.convert(sys.argv[1], 'fasta', sys.argv[1].replace('fasta', 'nexus')
              , 'nexus', alphabet=Alphabet.generic_dna)
