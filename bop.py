#!/usr/local/bin/python3

from glob import glob
from Bio import SeqIO
import argparse
import re

parser = argparse.ArgumentParser(description='''This program will add filename
into fasta sequence's id, and merge same gene files into single one.''')
parser.add_argument('path', default='./', help='target path, default is "./"')
arg = parser.parse_args()

filename = glob(arg.path+'*')
print(filename)
pattern = re.compile(r'.*/(.*)-(BOP\d{6}).fasta')
for i in filename:
    match =pattern.search(i)
    if match is not None:
        print(match.group(1), match.group(2))

