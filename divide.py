#!/usr/bin/python3

import sys
from Bio import SeqIO
from Bio import SearchIO
from Bio.Blast.Applications import NcbiblastnCommandline as nb
from os import makedirs
from os.path import exists

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

def main():
    """
    Usage:
    python3 divide.py fastq_file barcode_file primer_file
    This program may use large memory, be careful!
#    Before running this script, ensure you already put blast database file in
    current path. Also, it assumes that you have installed BLAST suite before. 
    To create the db file: 
    makeblastdb -in infile -out outfile -dbtype nucl

    Barcode file looks like this:
    ATACG BOP00001
    Primer file looks like this:
    rbcL rbcLF ATCGATCGATCGA
    From left to right, there are:
    1. gene name
    2. primer name
    3. primer sequence
    In this program, primer have common sequence whose length is 14, and
    barcode's is 5.
    """

    print(main.__doc__)
    if not exists('output'):
        makedirs('output')
    #read barcode file
    with open(sys.argv[2], 'r') as barcode_file:
        barcode_raw = barcode_file.read().split(sep='\n')
    barcode_raw.pop()
    barcode_list = [i.split() for i in barcode_raw]
    barcode = dict(barcode_list)
    #raise ValueError('test')
    #read primer file
    with open(sys.argv[3], 'r') as primer_file:
        primer_raw = primer_file.read().split(sep='\n')
    primer_raw.pop()
    primer = [i.split() for i in primer_raw]
    #large memory
    fastq_raw = list(SeqIO.parse(sys.argv[1], 'fastq'))
    divide_via_barcode(fastq_raw, barcode, primer)
    #convert fastq to fasta, then use BLAST to divide sequence via primer
    fasta_file = sys.argv[1].replace('fastq', 'fasta')
    #SeqIO.convert(sys.argv[1], 'fastq', fasta_file, 'fasta')

    fasta_raw = list(SeqIO.parse(fasta_file, 'fasta'))
    #blast(option)
    #parse(target)
    #output(target, option)


def divide_via_barcode(fastq_raw, barcode, primer):
    #change if necessary
    primer_adapter = 14
    barcode_length = 5
    total = len(fastq_raw)
    not_found = 0
    for record in fastq_raw:
        record_barcode = str(record.seq[:5])
        print(record_barcode)
        try:
            number = barcode[record_barcode]
        except:
            not_found += 1
        pass
    print(not_found, total)
    exit -1
if __name__ =='__main__':
    main()
