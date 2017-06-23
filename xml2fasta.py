#!/usr/bin/python3

from Bio import SearchIO
from Bio import SeqIO
from sys import argv


xml = SearchIO.parse(argv[1], 'blast-xml')
handle_tsv = open('{}.tsv'.format(argv[1]), 'w')
handle_tsv.write('Query id\tbitscore\tSpecies name\thit id\n')
for query in xml:
    if query.description != '':
        query.id = ''.join([query.id, query.description])
    if len(query) == 0:
        with open(argv[1]+'_not_found.log', 'a') as not_found:
            not_found.write('{} not found!\n'.format(query.id))
        continue
    handle = open('{}.fasta'.format(query.id), 'w')
    SeqIO.write(query[0][0].query, handle, 'fasta')
    for hit in query:
        for hsp in hit:
            species_name = hsp.hit.description.split(' ')
            if species_name[0].isupper():
                species_name = '_'.join(species_name[:3])
            else:
                species_name = '_'.join(species_name[:2])
            info = '{}\t{}\t{}\t{}{}\n'.format(
                query.id,
                hsp.bitscore,
                species_name,
                hsp.hit.id, hsp.hit.description)
            handle_tsv.write(info)
            hsp.hit.id = '{}|{}'.format(hsp.bitscore, hsp.hit.id)
            SeqIO.write(hsp.hit, handle, 'fasta')
