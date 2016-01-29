#!/usr/local/bin/python3

from Bio import SeqIO
from os import listdir
import argparse
import re

parser = argparse.ArgumentParser(description='''This program will add filename
into fasta sequence's id, and merge same gene files into single one.''')
parser.add_argument('--path', default='./', help='target path, default is "./"')
arg = parser.parse_args()

filename = listdir(arg.path)
pattern = re.compile(r'(.*)-(.*).fasta')
for i in filename:
    match =pattern.search(i)
    if match is not None:
        gene = match.group(1)
        sample = match.group(2)
        print(gene, sample)

