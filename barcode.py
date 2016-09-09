#!/usr/bin/python3

import argparse
import os
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline as nb
from glob import glob
from multiprocessing import cpu_count
from subprocess import run


def find_longest(fasta_files):
    avg_length = list()
    for fasta in fasta_files:
        length = list()
        raw = SeqIO.parse(fasta, 'fasta')
        for sequence in raw:
            length.append(len(sequence))
        avg_length.append([fasta, sum(length)/len(length)])
    avg_length.sort(key=lambda i:i[1])
    return  [i[0] for i in avg_length]


def get_sample(fasta_files, target):
    new_fasta_files = list()
    for fasta in fasta_files:
        output = fasta.replace('.fasta', '_{0}.fasta'.format(target))
        raw = SeqIO.parse(fasta, 'fasta')
        with open(output, 'w') as output_file:
            for n in range(target):
               SeqIO.write(next(raw), output_file, 'fasta')
        new_fasta_files.append(output)
    return new_fasta_files


def makeblastdb(db_file):
    db_name = db_file.replace('.fasta', '')
    run('makeblastdb -in {0} -title {1} -out {1} -dbtype nucl'.format(
        db_file, db_name), shell=True)
    return db_name


def blast(query_file, db_file, output_file='BLASTResult.xml'):
    """Here we use "max_hsps" to restrict only first hsp.
    """
    cmd = nb(num_threads=cpu_count(),
             query=query_file,
             db=db_file,
             task='blastn',
             outfmt=5,
             out=output_file)
    stdout, stderr = cmd()
    return output_file


def main():
    """This program will try to find out barcode to devide different species
    while ignore distinction among subspecies level."""
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('--path', default='.',
                        help='target path, default is present directory')
    parser.add_argument('--db', default=None, help='fasta file to make blast database, which contains longest sequence')
    parser.add_argument('--sample', default=None, type=int, help='sample numbers')
    parser.print_help()
    arg = parser.parse_args() 
    fasta_files = glob(arg.path+'*.fasta')
    if arg.sample is not None:
        fasta_files = get_sample(fasta_files, arg.sample)
    if arg.db is None:
        *query, db = find_longest(fasta_files)
    else:
        db = arg.db
        query = set(fasta_files) - db
    db_name = makeblastdb(db)
    blast_result = list()
    for fasta in query:
        result_file = fasta.replace('.fasta', '.xml')
        blast_result.append(blast(fasta, db_name, result_file))


if __name__ == '__main__':
    main()