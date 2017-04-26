#!/usr/bin/python3

from Bio import SeqIO
from sys import argv


handle = open(argv[1]+'.fasta', 'w')

for record in SeqIO.parse(argv[1], 'gb'):
    # order|family|organims(genus|species)
    organism = record.annotations['organism'].replace(' ', '_')
    order_family = record.annotations['taxonomy'][-3:-1]
    genus, *species  = organism.split('_')
    taxon = '{}|{}|{}|{}'.format(*order_family, genus, species)
    print(taxon)

