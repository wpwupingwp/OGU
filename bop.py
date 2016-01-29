#!/usr/local/bin/python3

from Bio import SeqIO
import argparse
import re
import os

parser = argparse.ArgumentParser(description='''This program will add filename
into fasta sequence's id, and merge same gene files into single one.''')
parser.add_argument('--path', default='./', help='target path, default is "./"')
arg = parser.parse_args()

filename = os.listdir(arg.path)
os.chdir(arg.path)
pattern = re.compile(r'(.*)-(.*).fasta')
gene = ''
sample = ''
for i in filename:
    match = pattern.search(i)
    if match is None:
        pass
    else:
        gene = match.group(1)
        handle = open('merge-'+gene, 'a')
        sample = match.group(2)
        fasta_file = SeqIO.parse(i, 'fasta')
        #print('{0} is not a legal fasta file.'.format(i))
        for sequence in fasta_file:
            sequence.id = re.sub(gene,'-'.join([gene, sample]), sequence.id)
            SeqIO.write(sequence, handle, 'fasta')
