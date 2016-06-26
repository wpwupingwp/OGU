#!/usr/local/bin/python3

from Bio import SeqIO
from time import process_time
import argparse
import re
import os

parser = argparse.ArgumentParser(description='''After process of divide.py, this program will add filename 
into fastq sequence's id, and merge same gene files into single one.''')
parser.add_argument('--path', default='./', help='target path, default is "./"')
arg = parser.parse_args()

filename = os.listdir(arg.path)
filename = [i for i in filename if i[-1] == 'q']
#only contain fastq files
os.chdir(arg.path)
pattern = re.compile(r'(.*)_(.*).fasta')
sample = ''
gene = ''
for i in filename:
    match = pattern.search(i)
    if match is None:
        pass
    else:
        gene = match.group(1)
        handle = open('merge-'+gene, 'a')
        sample = match.group(2)
        fasta_file = SeqIO.parse(i, 'fastq')
        for sequence in fasta_file:
            sequence.description = '-'.join([gene, sample, sequence.description])
#to avoid repeat replacement
            sequence.id = ''
            SeqIO.write(sequence, handle, 'fastq')
print('Finished in {:.3f}s.'.format(process_time()))
