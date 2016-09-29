#!/usr/bin/python3

import argparse
import os
from Bio import SearchIO, SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline as nb
from multiprocessing import cpu_count
from subprocess import run
from timeit import default_timer as timer


def divide_barcode(barcode_len, skip):
    # get barcode dict
    barcode = dict()
    with open(arg.barcode_file, 'r') as input_file:
        for line in input_file:
            if line.startswith('barcode'):
                continue
            line = line.split(sep=',')
            barcode[line[0]] = line[1].strip()
    # analyze input files
    divided_files = list()
    SEARCH_LEN = 20
    half = barcode_len // 2
    statistics = {'mismatch': 0, 'mode_wrong': 0, 'total': 0}
    fastq_raw = SeqIO.parse(arg.input, 'fastq')
    handle_miss = open(os.path.join(arg.output, 'barcode_miss.fastq'), 'w')
    head_file = os.path.join(arg.output, 'head.fasta')
    handle_fasta = open(head_file, 'w')
    for record in fastq_raw:
        statistics['total'] += 1
        # ignore wrong barcode
        if str(record.seq[:half]) != str(record.seq[half:barcode_len]):
            statistics['mode_wrong'] += 1
            continue
        record_barcode = [str(record.seq[:barcode_len]),
                          str(record.seq[:-(barcode_len + 1):-1])]
        if arg.strict:
            condition = (record_barcode[0] in barcode and
                         record_barcode[1] in barcode)
        else:
            condition = (record_barcode[0] in barcode)
        if condition:
            name = barcode[record_barcode[0]]
            output_file = os.path.join(arg.output, name)
            divided_files.append(output_file)
            with open(output_file, 'a') as handle:
                SeqIO.write(record, handle, 'fastq')
            handle_fasta.write('>{0}\n{1}\n'.format(
                record.description,
                record.seq[skip:skip + SEARCH_LEN]))
        else:
            SeqIO.write(record, handle_miss, 'fastq')
            statistics['mismatch'] += 1
    handle_miss.close()
    handle_fasta.close()
    return statistics, head_file, divided_files


def blast_and_parse(query_file, db_file):
    # Use blastn-short for primers.
    blast_result_file = os.path.join(arg.output, 'BlastResult.xml')
    cmd = nb(
        num_threads=cpu_count(),
        query=query_file,
        db=db_file,
        task='blastn-short',
        max_target_seqs=1,
        max_hsps=1,
        evalue=arg.evalue,
        outfmt=5,
        out=blast_result_file
    )
    stdout, stderr = cmd()
    # parse
    blast_result = SearchIO.parse(blast_result_file, 'blast-xml')
    parse_result = dict()
    for record in blast_result:
        # skip empty blast result item
        if len(record) == 0:
            continue
        else:
            tophit = record[0]
        query_info = '{0} {1}'.format(
            tophit[0][0].query_id,
            tophit[0][0].query_description)
        hit_info = tophit[0][0].hit.id
        parse_result[query_info] = hit_info
    return parse_result


def divide_gene(head_file, divided_files):
    # generate primer db
    primer = list()
    with open(arg.primer_file, 'r') as input_file:
        for line in input_file:
            if line.startswith('gene'):
                continue
            line = line.split(sep=',')
            primer.append([i.strip() for i in line])
    # join primer pairs
    join_seq = 'N'*15
    gene_list = list()
    primer_file = os.path.join(arg.output, 'primer.fasta')
    handle = open(primer_file, 'w')
    for index in range(0, len(primer) - 1, 2):
        left = primer[index][2][arg.primer_adapter:]
        right = primer[index + 1][2][arg.primer_adapter:]
        short_primer = ''.join([left, join_seq, right])
        name = primer[index][0]
        gene_list.append(name)
        handle.write('>{0}\n{1}\n'.format(name, short_primer))
    handle.close()
    # blast and parse
    db_name = primer_file.replace('.fasta', '')
    run('makeblastdb -in {0} -out {1} -dbtype nucl'.format(
        primer_file, db_name), shell=True)
    blast_result = blast_and_parse(head_file, db_name)
    # split and count
    sample_count = {i: 0 for i in divided_files}
    gene_count = {i: 0 for i in gene_list}
    for fastq_file in divided_files:
        records = SeqIO.parse(fastq_file, 'fastq')
        for record in records:
            gene = record.description
            if gene in blast_result:
                sample_count[fastq_file] += 1
                gene_count[blast_result[gene]] += 1
                barcode = os.path.basename(fastq_file).replace('.fastq', '')
                record.id = '|'.join([barcode, record.id])
                handle = open('{0}_{1}.fastq'.format(
                    fastq_file.replace('.fastq', ''), blast_result[gene]), 'a')
                SeqIO.write(record, handle, 'fastq')
                if not arg.no_merge_gene:
                    handle_gene = open(os.path.join(
                        arg.output, gene+'.fastq'), 'a')
                    SeqIO.write(record, handle_gene, 'fastq')
    return sample_count, gene_count


def print_stats(barcode, sample, gene):
    with open(os.path.join(arg.output, 'barcode_info.csv'), 'w') as handle:
        for record in barcode.items():
            handle.write('{0},{1} \n'.format(*record))
    with open(os.path.join(arg.output, 'sample_info.csv'), 'w') as handle:
        for record in sample.items():
            handle.write('{0},{1} \n'.format(*record))
    with open(os.path.join(arg.output, 'gene_info.csv'), 'w') as handle:
        for record in gene.items():
            handle.write('{0},{1} \n'.format(*record))


def main():
    # todo
    # implement mode
    # check input
    """
    Step 1, divide data by barcode.
    Step 2, divide data by primer via BLAST.
    Ensure that you have installed BLAST suite before.
    Make sure you don't miss the first line.
    Barcode file looks like this:
    ATACG,BOP00001
    Primer file looks like this:
    gene,primer,sequence,direction
    rbcL,rbcLF,ATCGATCGATCGA
    rbcL,rbcLR,TACGTACGTACG
    Make sure you don't miss the first line and the capitalization.
    To get these two files, save your excel file as csv file.  Be carefull
    of the order of  each pair of primers.
    From left to right, there are:
    1. gene name
    2. primer name
    3. primer sequence
    4. forward/backward
    """
    print(main.__doc__)
    start_time = timer()
    parser = argparse.ArgumentParser()
    parser.add_argument('--barcode_length', default=10, type=int,
                        help='length of barcode')
    parser.add_argument('--primer_adapter', default=14, type=int,
                        help='length of primer_adapter, typical 14 for AFLP')
    parser.add_argument('-b', dest='barcode_file',
                        help='csv file containing barcode info')
    parser.add_argument('-p', dest='primer_file',
                        help='csv file containing primer info')
    parser.add_argument('-e', dest='evalue', default=1e-5, type=float,
                        help='evalue for BLAST')
    parser.add_argument('-s', '--strict', action='store_true',
                        help='''if set nostrict, it will only consider
                        barcode on the head; if not, consider head and tail''')
    parser.add_argument('-m', dest='mode', default='5-2',
                        help='''barcode mode, default value is 5-2, i.e.,
                        barcode with length 5 repeated twice''')
    parser.add_argument('--no_merge_gene', action='store_true',
                        help='merge output files by gene')
    parser.add_argument('input', help='input file, fastq format')
    parser.add_argument('-o', dest='output', default='out', help='output path')
    global arg
    arg = parser.parse_args()
    skip = arg.barcode_length + arg.primer_adapter
    if not os.path.exists(arg.output):
        os.mkdir(arg.output)

    barcode_info, head_file, divided_files = divide_barcode(
        arg.barcode_length, skip)
    sample_count, gene_count = divide_gene(head_file, divided_files)
    print_stats(barcode_info, sample_count, gene_count)

    end_time = timer()
    print('Finished with {0:.3f}s. You can find results in {1}.\n'.format(
        end_time-start_time, arg.output))

if __name__ == '__main__':
    main()
