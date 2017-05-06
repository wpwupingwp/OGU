#!/usr/bin/python3

from glob import glob
from Bio import SeqIO

info = dict()
with open('./gb_info.csv', 'r') as raw:
    for raw_line in raw:
        line = raw_line.strip()
        line = line.split(sep=',')
        info[line[0]] = line

P = 'v3/*.fasta'
file_list = glob(P)
for i in file_list:
    with open(i+'.rename', 'w') as output:
        for record in SeqIO.parse(i, 'fasta'):
            old_id = record.description
            accession, gene = old_id.split('_')
            try:
                accession_id, *taxon = info[accession]
                record.description = '|'.join([gene, *taxon, accession_id])
            except:
                print(record.description)
                record.description = '|'.join([gene, accession])
            record.id = ''
            SeqIO.write(record, output, 'fasta')
