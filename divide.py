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


def step1():
    """
    Divide raw data via barcode.
    ID were added into the beginning of sequence id."""
    print(step1.__doc__)
    barcode = get_barcode_dict()
    fastq_raw = SeqIO.parse(sys.argv[1], 'fastq')
    total = 0
    not_found = 0
    handle = open('step1.fastq', 'w')
    handle2 = open('step1_miss.fastq', 'w')
    for record in fastq_raw:
        total += 1
        record_barcode = str(record.seq[:5])
        try:
            new_id = barcode[record_barcode]
        except:
            SeqIO.write(record, handle2, 'fastq')
            not_found += 1
            continue
        new_record = SeqRecord(
            id='|'.join([new_id, record.id]),
            description=record.description,
            seq=record.seq,
            letter_annotations=record.letter_annotations
        )
        SeqIO.write(new_record, handle, 'fastq')
    SeqIO.convert('step1.fastq', 'fastq', 
                  'step1.fasta', 'fasta')

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


def write_fasta(primer_list, skip):
    handle = open('primer.fasta', 'w')
    join_seq = 'NNNNNNNNNNNNNNN'
    for index in range(0, len(primer_list)-1, 2):
        left = primer_list[index][2][skip:]
        right = primer_list[index+1][2][skip:]
        short_primer = ''.join([left, join_seq, right])
        name = str(primer_list[index][0],)
        handle.write(''.join([
            '>', 
            name, 
            '\n',
            short_primer,
            '\n'
        ]))
    handle.close()


def step2():
    """
    BLAST fastq in first step against primer database.  Next, use BLAST 
    result to divide again."""
    print(step2.__doc__)
    barcode_length = 5
    primer_adapter = 14
    skip = barcode_length + primer_adapter
    primer_list = get_primer_list()
    write_fasta(primer_list, skip)
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
    if not exists('out'):
        makedirs('out')
    miss_step1, total = step1()
    print('''
    Step1 results:
    Total: {0} reads
    unrecognize {1} reads 
    {2:3f} percent'''.format(total, miss_step1, miss_step1/total))
    step2()

if __name__ =='__main__':
    main()
