#!/usr/bin/python3

from Bio import SeqIO
from Bio import SearchIO
from Bio.Blast.Applications import NcbiblastnCommandline as nb
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from glob import glob
from multiprocessing import cpu_count
from os import makedirs
from os.path import exists
from subprocess import call
import sys


def parse(target):
    blast_result = list(SearchIO.parse('BlastResult.xml', 'blast-xml'))
    for record in _result:
        if len(record) == 0:
            continue
        else:
            tophit = record[0]
        target.append([tophit[0][0].query, tophit[0][0].hit])


def get_barcode_dict():
    with open(sys.argv[2], 'r') as input_file:
        barcode_raw = input_file.read().split(sep='\n')
    barcode_raw.pop(0)
    barcode_raw.pop(-1)
    barcode_list = [i.split(sep=',') for i in barcode_raw]
    barcode_dict = dict(barcode_list)
    print(barcode_dict)
    return barcode_dict


def step1(skip):
    """
    Divide raw data via barcode."""
    print(step1.__doc__)
    barcode = get_barcode_dict()
    fastq_raw = SeqIO.parse(sys.argv[1], 'fastq')
    total = 0
    not_found = 0
    handle_miss = open('step1_miss.fastq', 'w')
    handle_fasta = open('step1.fasta', 'w')
    for record in fastq_raw:
        total += 1
        record_barcode = str(record.seq[:5])
        try:
            name = barcode[record_barcode]
        except:
            SeqIO.write(record, handle_miss, 'fastq')
            not_found += 1
            continue
        handle = open(name+'.fastq', 'a')
        SeqIO.write(record, handle, 'fastq')
    #only use head to blast
        handle_fasta.write(''.join([
            '>', name, '|', record.description, '\n',
            str(record.seq[skip:skip+20]), '\n'
        ]))
    handle.close()
    handle_miss.close()
    handle_fasta.close()
    return (not_found, total)

def blast(query_file, database):
    cmd = nb(
        num_threads=cpu_count(),
        query=query_file,
        db=database,
        task='blastn-short', 
        evalue=1e-5, 
        outfmt=5, 
        out='out/BlastResult.xml'
    )
    stdout, stderr = cmd()

def get_primer_list():
    with open(sys.argv[3], 'r') as input_file:
        primer_raw = input_file.read().split(sep='\n')
    primer_raw.pop(0)
    primer_raw.pop(-1)
    primer_list = [i.split(sep=',') for i in primer_raw]
    return primer_list


def write_fasta(primer_list, primer_adapter):
    handle = open('primer.fasta', 'w')
    for index in range(len(primer_list)):
        short_primer = primer_list[index][2][primer_adapter:]
        name = ''.join([primer_list[index][0], '-', str(index)])
        handle.write(''.join([
            '>', name, '\n',
            short_primer, '\n'
        ]))
    handle.close()

def write_fasta_2(primer_list, primer_adapter):
    handle = open('primer.fasta', 'w')
    join_seq = 'NNNNNNNNNNNNNNN'
    for index in range(0, len(primer_list)-1, 2):
        left = primer_list[index][2][primer_adapter:]
        right = primer_list[index+1][2][primer_adapter:]
        short_primer = ''.join([left, join_seq, right])
        name = primer_list[index][0]
        handle.write(''.join([
            '>', name, '\n',
            short_primer, '\n'
        ]))
    handle.close()


def step2(primer_adapter):
    """
    BLAST fastq in first step against primer database.  Next, use BLAST 
    result to divide again."""
    print(step2.__doc__)
    primer_list = get_primer_list()
    write_fasta(primer_list, primer_adapter)
    call('makeblastdb -in primer.fasta -out primer -dbtype nucl',
         shell=True) 
    blast('step1.fasta', 'primer')


def main():
    """
    Usage:
    python3 divide.py fastqFile barcodeFile primerFile
    Step 1, divide data by barcode. Step 2, divide data by primer via BLAST.
    Ensure that you have installed BLAST suite before. 
    Barcode file looks like this:
    ATACG,BOP00001
    Primer file looks like this:
    rbcL,rbcLF,ATCGATCGATCGA
    rbcL,rbcLR,TACGTACGTACG
    To get these two files, save your excel file as csv file.  Be carefull 
    of the order of  each pair of primers.
    From left to right, there are:
    1. gene name
    2. primer name
    3. primer sequence
    In this program, primer have common sequence whose length is 14, and
    barcode's is 5. Change them if necessary(in step2)."""
    print(main.__doc__)
    barcode_length = 5
    primer_adapter = 14
    skip = barcode_length + primer_adapter
    if not exists('out'):
        makedirs('out')
    miss_step1, total = step1(skip)
    print('''
    Step1 results:
    Total: {0} reads
    unrecognize {1} reads 
    {2:3f} percent'''.format(total, miss_step1, miss_step1/total))
    step2(primer_adapter)

if __name__ =='__main__':
    main()
