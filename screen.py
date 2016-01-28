from Bio import SeqIO
import re
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename', help='fasta file')
parser.add_argument('-l', '--length', type=int, dest='length', default=50, help='minimal sequence length')
parser.add_argument('-c', '--cover', type=int, dest='cover', default=1, help='minimal coverage')
arg = parser.parse_args()
print(arg.length, arg.cover)