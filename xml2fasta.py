#!/usr/bin/python3

from Bio import SearchIO
from Bio import SeqIO
from sys import argv


xml = SearchIO.parse(argv[1], 'blast-xml')
for query in xml:
    if len(query) == 0:
        print(query.id)
    for hit in query:
        print(len(hit))
        for hsp in hit:
            print(hsp.bitscore)
