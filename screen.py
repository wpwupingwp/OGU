#!/usr/bin/python3

# Screen contigs according to coverage and length
import re
import argparse
from Bio import SeqIO
from time import process_time

parser = argparse.ArgumentParser()
parser.add_argument('filename', help='fasta file')
parser.add_argument('-l', '--length', type=int, dest='length',
                    default=0, help='minimal sequence length')
parser.add_argument('-c', '--cover', type=float, dest='cover',
                    default=0, help='minimal coverage')
arg = parser.parse_args()

with open(arg.filename, 'r') as input_file:
    raw = list(SeqIO.parse(input_file, 'fasta'))
output_file = arg.filename + '.filtered'
handle = open(output_file, 'a')

pattern = re.compile(r'length_(\d+)_cov_(\d+.\d+)')
count = 0
for i in raw:
    m = pattern.search(i.id)
    if m is not None:
        count += 1
        length = int(m.group(1))
        cover = float(m.group(2))
        if length > arg.length and cover > arg.cover:
            SeqIO.write(i, handle, 'fasta')
            print('No.{0}\tlength: {1}\tCover: {2}'.format(
                count, length, cover))
print('Finished within {:.3f}s.'.format(process_time()))
