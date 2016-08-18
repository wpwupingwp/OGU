#!/usr/bin/python3

import argparse
from Bio import SeqIO

parameters = argparse.ArgumentParser(
    description='Generate fasta file with specific id format.')
parameters.add_argument('-i', '--input', help='gb format file')
parameters.add_argument('-o', '--output', help='outpu file')
parameters.print_help()
arg = parameters.parse_args()
for sequence in SeqIO.parse(arg.input, 'genbank'):
    for feature in sequece.features:
        pass
