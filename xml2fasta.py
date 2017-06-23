#!/usr/bin/python3

from Bio import SearchIO
from Bio import SeqIO
import argparse


def main():
    arg = argparse.ArgumentParser()
    arg.add_argument('input', help='input BLAST result (xml format)')
    arg.add_argument('-s', '--simple', action='store_true',
                     help='only handle first hsp')
    arg.add_argument('-ss', '--very_simple', action='store_true',
                     help='only handle first hit')
    arg = arg.parse_args()

    xml = SearchIO.parse(arg.input, 'blast-xml')
    handle_tsv = open('{}.tsv'.format(arg.input), 'w')
    handle_tsv.write('Query id\tbitscore\tSpecies name\thit id\n')
    for query in xml:
        if query.description != '':
            query.id = ''.join([query.id, query.description])
        if len(query) == 0:
            with open(arg.input+'_not_found.log', 'a') as not_found:
                not_found.write('{} not found!\n'.format(query.id))
            continue
        handle = open('{}.fasta'.format(query.id), 'w')
        SeqIO.write(query[0][0].query, handle, 'fasta')
        for hit in query:
            for hsp in hit:
                species_name = hsp.hit.description.split(' ')
                if species_name[0].isupper():
                    species_name = '{}_{}_{}'.format(
                        *species_name[1:3], species_name[0].replace(':', ''))
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
                if arg.simple or arg.very_simple:
                    break
            if arg.very_simple:
                break


if __name__ == '__main__':
    main()
