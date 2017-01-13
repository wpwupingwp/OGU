#!/usr/bin/python3

from Bio import SeqIO
from Bio import SearchIO
from Bio.Blast.Applications import NcbiblastnCommandline as nb
from functools import wraps
from timeit import default_timer as timer
from multiprocessing import cpu_count
from os import path, mkdir
from subprocess import call
import argparse


def print_time(function):
    @wraps(function)
    def wrapper(*args, **kargs):
        start = timer()
        result = function(*args, **kargs)
        end = timer()
        print('The function {0} costed {1:.3f}s.'.format(
            function.__name__, end-start))
        return result
    return wrapper


@print_time
def blast(query_file, db_file, output_file='BLASTResult.xml'):
    """Here we use "max_hsps" to restrict only first hsp.
    """
    call('makeblastdb -in {0} -out {0} -dbtype nucl'.format(db_file),
         shell=True)
    cmd = nb(num_threads=cpu_count(),
             query=query_file,
             db=db_file,
             task='blastn',
             evalue=arg.evalue,
             max_hsps=1,
             max_target_seqs=1,
             outfmt=5,
             out=output_file)
    stdout, stderr = cmd()
    return output_file


@print_time
def parse(blast_result):
    raw = list()
    result = SearchIO.parse(blast_result, 'blast-xml')
    for query in result:
        for hit in query:
            for hsp in hit:
                cover = (hsp.query_end - hsp.query_start) / (hsp.seq.__len__)
                if cover < arg.cover:
                    continue
                line = [hsp.query, hsp.hit]
                raw.append(line)
    return raw


def output():
    SeqIO.write()


def main():
    parameters = argparse.ArgumentParser(
        description='Filter contigs according to given genome sequence.')
    parameters.add_argument('-i', '--input', help='fasta format genome file')
    parameters.add_argument('-r', '--reference',
                            help='reference genome(fasta format)')
    parameters.add_argument('-c', '--cover', type=float, default=0.5,
                            help='coverage of blast')
    parameters.add_argument('-e', '--evalue', type=float, default=1e-5,
                            help='evalue for BLAST')
    parameters.add_argument('-o', '--output', default='out',
                            help='output directory')
    parameters.print_help()
    global arg
    arg = parameters.parse_args()
    if not path.exists(arg.output):
        mkdir(arg.output)


if __name__ == '__main__':
    main()
