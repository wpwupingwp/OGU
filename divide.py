#!/usr/bin/python3

from Bio import SeqIO
from Bio import SearchIO
from Bio.Blast.Applications import NcbiblastnCommandline as nb
from Bio.SeqRecord import SeqRecord
from glob import glob
from os import makedirs
from os.path import exists
from subprocess import call
from multiprocessing import cpu_count
import sys


def parse(target):
    blast_result = list(SearchIO.parse('BlastResult.xml', 'blast-xml'))
    for record in blast_result:
        if len(record) == 0:
            continue
        else:
            tophit = record[0]
        target.append([tophit[0][0].query, tophit[0][0].hit])


def step1():
    """
    Divide raw data via barcode.
    Results were placed in 'output/'."""
    print(step1.__doc__)

    with open(sys.argv[2], 'r') as barcode_file:
        barcode_raw = barcode_file.read().split(sep='\n')
    barcode_raw.pop()
    barcode_list = [i.split() for i in barcode_raw]
    barcode = dict(barcode_list)
    fastq_raw = SeqIO.parse(sys.argv[1], 'fastq')
    total = 0
    not_found = 0
    found = list()
    for record in fastq_raw:
        total += 1
        record_barcode = str(record.seq[:5])
        try:
            number = barcode[record_barcode]
            output_file = ''.join(['output/', number, '.fastq'])
            handle = open(output_file, 'a')
            SeqIO.write(record, handle, 'fastq')
            found.append(record)
        except:
            not_found += 1
    exit -1
    return (not_found, total)

def blast():
    cmd = nb(
    #change if necessary
        num_threads=4,
        query=query_file,
        db=sys.argv[2], 
        task='blastn-short', 
        evalue=0.001, 
        outfmt=5, 
        # xml format
        out='BlastResult.xml'
    )
    stdout, stderr = cmd()
    pass


def step2():
    """
    At first, create blast database file with makeblastdb.
    Database will be put in 'output/', and named 'primer'.
    Also convert fastq file to fasta file.
    Then BLAST fastq in first step against primer database.
    Next, use BLAST result to divide again."""
    print(create_blast_db.__doc__)
    barcode_length = 5
    primer_adapter = 14
    skip = barcode_length + primer_adapter

    with open(sys.argv[3], 'r') as primer_file:
        primer_raw = primer_file.read().split(sep='\n')
    primer_raw.pop(0)
    primer_raw.pop(-1)
    primer_list = [i.split() for i in primer_raw]

    primer_fasta = list()
    for index in range(0, len(primer_list), 2):
        left = primer[index][2][skip:]
        right = primer[index+1][2][skip:]
        short_primer = 'NNNNN'.join([left, right])
        sequence = SeqRecord(
            id=primer[index][0],
            description='',
            seq=short_primer
        )
        primer_fasta.append(sequence)
    SeqIO.write(primer_fasta, 'output/primer.fasta', 'fasta')

    fasta_file = sys.argv[1].replace('fastq', 'fasta')
    SeqIO.convert(sys.argv[1], 'fastq', fasta_file, 'fasta')
    call('makeblastdb -in output/primer.fasta -out output/primer -dbtype nucl') 


def main():
    """
    Usage:
    python3 divide.py fastq_file barcode_file primer_file
    Step 1, divide data by barcode.
    Step 2, divide data by primer via BLAST.
    Ensure that you have installed BLAST suite before. 
    Barcode file looks like this:
    ATACG BOP00001
    Primer file looks like this:
    rbcL rbcLF ATCGATCGATCGA
    rbcL rbcLR TACGTACGTACG
    Be carefull of the order of  each pair of primers.
    From left to right, there are:
    1. gene name
    2. primer name
    3. primer sequence
    In this program, primer have common sequence whose length is 14, and
    barcode's is 5. Change them if necessary(in step2_prepare).  """
    print(main.__doc__)
    if not exists('output'):
        makedirs('output')
    miss_step1, total = step1()
    step2()

if __name__ =='__main__':
    main()
