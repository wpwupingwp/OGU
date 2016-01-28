from Bio import SeqIO
import re
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename', help='fasta file')
parser.add_argument('-l', '--length', type=int, dest='length', default=0, help='minimal sequence length')
parser.add_argument('-c', '--cover', type=float, dest='cover', default=0, help='minimal coverage')
arg = parser.parse_args()

with open(arg.filename, 'r') as input_file:
    raw = list(SeqIO.parse(input_file, 'fasta'))
output_file = arg.filename + '.filtered'
handle = open(output_file, 'w')
for i in raw:
    m = re.search(r'length_(\d+)_cov_(\d+.\d)_', i.id)
    if m is not None:
        length = int(m.group(1))
        cover = float(m.group(2))
        if length > arg.length and cover > arg.cover:
            SeqIO.write(i, handle, 'fasta')
            print(length,cover)
