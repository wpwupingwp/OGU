#!/usr/bin/python3

from Bio import SearchIO
from Bio import SeqIO
from sys import argv


xml = SearchIO.parse(argv[1], 'blast-xml')
for query in xml:
    if len(query) == 0:
        print('{} not found!'.format(query.id))
        continue
    print(query.id)
    handle = open('{}_{}.fasta'.format(query.id, query.description), 'w')
    SeqIO.write(query[0][0].query, handle, 'fasta')
    for hit in query:
        for hsp in hit:
            print(hsp.bitscore, hsp.hit.id, hsp.hit.description)
            hsp.hit.id = '{}|{}'.format(hsp.bitscore, hsp.hit.id)
            SeqIO.write(hsp.hit, handle, 'fasta')
