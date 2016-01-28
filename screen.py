from Bio import SeqIO
import re
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename', help='fasta file')
parser.add_argument('-l', '--length', type=int, dest='length', default=0, help='minimal sequence length')
parser.add_argument('-c', '--cover', type=float, dest='cover', default=0, help='minimal coverage')
arg = parser.parse_args()

raw = list(SeqIO.parse(arg.filename, 'fasta'))
for i in raw:
    m = re.search(r'length_(\d+)_cov_(\d+.\d)_', i.id)
    if m is not None:
        length = m.group(1)
        cover = m.group(2)
        print(length,cover)
