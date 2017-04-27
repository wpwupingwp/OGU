#!/usr/bin/python3

from glob import glob
from Bio import SeqIO


info = dict()
with open('info.csv', 'r') as raw:
    for line in raw:
        line = line.strip().split(',')
        key = line[-1]
        info[key] = line

file_list = glob('*.fasta')
for fasta in file_list:
    handle = open(fasta+'.rename', 'w')
    for seq in SeqIO.parse(fasta, 'fasta'):
        # >A11|ycf1|NODE
        library, gene, *_ = seq.id.split('|')
        try:
            seq.id = '|'.join([gene, *info[library]])
        except:
            print(library)
        seq.description = ''
        SeqIO.write(seq, handle, 'fasta')
        # print(seq.id)
        # print(seq.description)
