#!/usr/bin/python3

import sys
from Bio import SeqIO
from Bio import SearchIO
from Bio.Blast.Applications import NcbiblastnCommandline as nb
from os import makedirs
from os.path import exists

def blast(option):
    if option == '1':
        query_file = './output/all.fasta'
    else:
        query_file = sys.argv[1]
    cmd = nb(
    #  num_threads=8,
        query=query_file,
        db=sys.argv[2], 
        task='megablast', 
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
    """\n
    Usage:
    python3 divide.py fastq_file barcode_file primer_file
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
    3. primer sequence"""

    print(main.__doc__)
    if not exists('output'):
        makedirs('output')
    #read barcode file
    with open(sys.argv[2], 'r') as barcode_file:
        barcode_raw = barcode_file.read().split(sep='\n')
    barcode_raw.pop()
    barcode = [i.split() for i in barcode_raw]
    #read primer file
    with open(sys.argv[3], 'r') as primer_file:
        primer_raw = primer_file.read().split(sep='\n')
    primer_raw.pop()
    primer = [i.split() for i in primer_raw]
    print(primer)
    raise ValueError('test')

    fasta_file = sys.argv[1].replace('fastq', 'fasta')
    SeqIO.convert(sys.argv[1], 'fastq', fasta_file, 'fasta')
    blast(option)
    parse(target)
    output(target, option)

if __name__ =='__main__':
    main()
