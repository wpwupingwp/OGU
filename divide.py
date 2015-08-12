#!/usr/bin/python3

import sys
from Bio import SeqIO
from Bio import SearchIO
from Bio.Blast.Applications import NcbiblastnCommandline as nb
from Bio.SeqRecord import SeqRecord
from os import makedirs
from os.path import exists
from subprocess import call

def blast(option):
    cmd = nb(
    #  num_threads=8,
        query=query_file,
        db=sys.argv[2], 
        task='blastn-short', 
        evalue=0.001, 
        outfmt=5, 
        # xml format
        out='BlastResult.xml'
    )
    stdout, stderr = cmd()
    return 

def parse(target):
    blast_result = list(SearchIO.parse('BlastResult.xml', 'blast-xml'))
    for record in blast_result:
        if len(record) == 0:
            continue
        else:
            tophit = record[0]
        target.append([tophit[0][0].query, tophit[0][0].hit])

def output(target, option):
    if option == '1':
        for record in target:
            gene = record[0].id.split(sep='|')[-1]
            output_file = ''.join([
                'output/',
                sys.argv[2], 
                '-',
                gene, 
                '.fasta'
            ])
            rename_seq = SeqRecord(
                seq=record[1].seq, 
                id='|'.join([
                    gene,
                    sys.argv[1],
                    record[1].id
                ]),
                description=''
            )
            SeqIO.write(rename_seq, output_file, 'fasta')
    else:
        output_file = open('output/' + sys.argv[1] + '-filtered.fasta', 'w' )
        contig_id = {i[0].id for i in target} 
        query_file = SeqIO.parse(sys.argv[1], 'fasta')
        for record in query_file:
            if record.id in contig_id:
                SeqIO.write(record, output_file, 'fasta')


def step1():
    """Divide raw data via barcode."""
    print(step1.__doc__)

    with open(sys.argv[2], 'r') as barcode_file:
        barcode_raw = barcode_file.read().split(sep='\n')
    barcode_raw.pop()
    barcode_list = [i.split() for i in barcode_raw]
    barcode = dict(barcode_list)
    fastq_raw = SeqIO.parse(sys.argv[1], 'fastq')
    total = 0
    not_found = 0
    for record in fastq_raw:
        total += 1
        record_barcode = str(record.seq[:5])
        try:
            number = barcode[record_barcode]
            output_file = ''.join(['output/', number, '.fastq'])
            handle = open(output_file, 'a')
            SeqIO.write(record, handle, 'fastq')
        except:
            not_found += 1
    return (not_found, total)

def create_blast_db(barcode_length, primer_adapter)

    #read primer file
    skip = barcode_length + primer_adapter
    with open(sys.argv[3], 'r') as primer_file:
        primer_raw = primer_file.read().split(sep='\n')
    primer_raw.pop()
    primer = [i.split() for i in primer_raw]
    primer_fasta = list()
    for index in range(0, len(primer), 2):
        primer_short_left = primer[index][2][skip:]
        primer_short_right = primer[index+1][2][skip:]
        primer_short = 'NNNNN'.join([primer_short_left, primer_short_right])
        primer_fasta.append([
            id=primer[index][0],
            description='',
            seq=primer_short
        ])


    fasta_file = sys.argv[1].replace('fastq', 'fasta')
    SeqIO.convert(sys.argv[1], 'fastq', fasta_file, 'fasta')
    #create blast database file
    call(['makeblastdb -in', sys.argv[3], '-out primer -dbtype nucl']) 


def divide_via_primer(barcode_length, primer_adapter):
def main():
    """
    Usage:
    python3 divide.py fastq_file barcode_file primer_file
    Ensure that you have installed BLAST suite before. 
    Barcode file looks like this:
    ATACG BOP00001
    Primer file looks like this:
    rbcL rbcLF ATCGATCGATCGA
    From left to right, there are:
    1. gene name
    2. primer name
    3. primer sequence
    In this program, primer have common sequence whose length is 14, and
    barcode's is 5. Change them if necessary
    """
    barcode_length = 5
    primer_adapter = 14
    print(main.__doc__)
    if not exists('output'):
        makedirs('output')

    miss_step1, total = step1()
    create_blast_db(barcode_length, primer_adapter)
    step2(barcode_length, primer_adapter)
    #blast(option)
    #parse(target)
    #output(target, option)


if __name__ =='__main__':
    main()
