#!/usr/bin/python3

from Bio import SeqIO
from Bio import SearchIO
from Bio.Blast.Applications import NcbiblastnCommandline as nb
from multiprocessing import cpu_count
from os import makedirs
from os.path import exists
from subprocess import call
import sys
from glob import glob


def get_barcode_dict():
    with open(sys.argv[2], 'r') as input_file:
        barcode_raw = input_file.read().split(sep='\n')
    barcode_raw.pop(0)
    barcode_raw.pop(-1)
    barcode_list = [i.split(sep=',') for i in barcode_raw]
    barcode_dict = dict(barcode_list)
    return barcode_dict


def step1(skip):
    """
    Divide raw data via barcode.
    In this case, it searches primers in first 20bp sequence of reads, you may
    edit it."""
    print(step1.__doc__)
    search_len = 20
    barcode = get_barcode_dict()
    fastq_raw = SeqIO.parse(sys.argv[1], 'fastq')
    total = 0
    not_found = 0
    handle_miss = open('step1_miss.fastq', 'w')
    handle_fasta = open('step1.fasta', 'w')
    for record in fastq_raw:
        total += 1
        record_barcode = [str(record.seq[:5]), str(record.seq[-5::-1])]
        if record_barcode[0] in barcode and record_barcode[1] in barcode:
            name = barcode[record_barcode[0]]
            output_file = 'out/{0}'.format(name)
            with open(output_file, 'a') as handle:
                SeqIO.write(record, handle, 'fastq')
            handle_fasta.write('>{0}\n{1}\n'.format(
                record.description, 
                record.seq[skip:skip + search_len]))
        #elif record_barcode[1] in barcode:
        #    name = barcode[record_barcode[1]]
        #    output_file = 'out/'.format(name)
        #    with open(output_file, 'a') as handle:
        #        SeqIO.write(record, handle, 'fastq')
        #    handle_fasta.write('>{0}\n{1}\n'.format(
        #        record.description,
        #        record.seq[-(skip + search_len)::-1]))
        else:
            SeqIO.write(record, handle_miss, 'fastq')
            not_found += 1
    handle_miss.close()
    handle_fasta.close()
    return not_found, total


def get_primer_list():
    with open(sys.argv[3], 'r') as input_file:
        primer_raw = input_file.read().split(sep='\n')
    primer_raw.pop(0)
    primer_raw.pop(-1)
    primer_list = [i.split(sep=',') for i in primer_raw]
    return primer_list


def write_fasta(primer_list, primer_adapter):
    handle = open('primer.fasta', 'w')
    join_seq = 'NNNNNNNNNNNNNNN'
    gene_list = list()
    for index in range(0, len(primer_list) - 1, 2):
        left = primer_list[index][2][primer_adapter:]
        right = primer_list[index + 1][2][primer_adapter:]
        short_primer = ''.join([left, join_seq, right])
        name = primer_list[index][0]
        gene_list.append(name)
        handle.write('>{0}\n{1}\n'.format(name,short_primer))
    handle.close()
    return gene_list


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


def parse_blast():
    parse_result = list()
    blast_result = SearchIO.parse('out/BlastResult.xml', 'blast-xml')
    for record in blast_result:
        if len(record) == 0:
            continue
        else:
            tophit = record[0]
        query_info = '{0} {1}'.format(
            tophit[0][0].query_id,
            tophit[0][0].query_description)
        hit_info = tophit[0][0].hit.id
        parse_result.append([query_info, hit_info])
    parse_result = dict(parse_result)
    return parse_result


def step2(primer_adapter):
    """
    Step 2:
    BLAST fastq in first step against primer database."""
    print(step2.__doc__)
    primer_list = get_primer_list()
    gene_list = write_fasta(primer_list, primer_adapter)
    call('makeblastdb -in primer.fasta -out primer -dbtype nucl',
         shell=True)
    blast('step1.fasta', 'primer')
    blast_result = parse_blast()
    return blast_result, gene_list


def step3(blast_result, file_list, gene_list):
    """
    Step 3:
    First, according BLAST result, split fastq files generated in step1, then
    assembly."""
    print(step3.__doc__)
    count_sample = {i: 0 for i in file_list}
    count_gene = {i: 0 for i in gene_list}
    for fastq_file in file_list:
        records = SeqIO.parse(fastq_file, 'fastq')
        for record in records:
            gene = record.description
            if gene in blast_result:
                count_sample[fastq_file] += 1
                count_gene[blast_result[gene]] += 1
                handle = open(
                    '{0}_{1}'.format(fastq_file, blast_result[gene]), 
                    'a')
                SeqIO.write(record, handle, 'fastq')
    return count_sample, count_gene


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
    unrecognised {1} reads
    {2:3f} percent'''.format(total, miss_step1, miss_step1 / total))
    blast_result, gene_list = step2(primer_adapter)
    file_list = glob('out/B*')
    count_sample, count_gene = step3(blast_result, file_list, gene_list)
    count_sample = list(count_sample.items())
    count_gene = list(count_gene.items())
    with open('count_sample', 'w') as handle:
        for i in count_sample:
            handle.write('{0} {1} \n'.format(i[0], i[1]))
    with open('count_gene', 'w') as handle:
        for i in count_gene:
            handle.write('{0} {1} \n'.format(i[0], i[1]))
            # wrong


if __name__ == '__main__':
    main()
