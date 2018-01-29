#!/usr/bin/python3

from Bio import SeqIO
from sys import argv


def generate_fastq():
    # for Sequencher contig output
    seq_file, qual_file, fastq_file = argv[1:4]
    with open(seq_file, 'r') as seq, open(qual_file, 'r') as qual:
        fastq = SeqIO.QualityIO.PairedFastaQualIterator(seq, qual)
        SeqIO.write(fastq, fastq_file, 'fastq')
    return fastq_file


def get_id_info(record):
    raw_id = record.id
    print(raw_id)
    name2, tm, sample, cov, sum_bitscore_raw, location = raw_id.split('-')
    cov = cov.strip('%')
    return [name2, tm, sample, float(cov)/100, sum_bitscore_raw,
            int(float(location))]


def get_resolution(start, end):
    with open(argv[4], 'r') as resolution_file:
        raw = resolution_file.readlines()
        resolution = [i.split('\t')[1] for i in raw]
        if start > end:
            start, end = end, start
        fragment = resolution[start]
    return float(fragment)/100


def main():
    print('Usage:')
    print('python3 get_primer_info.py SeqFile QualFile'
          ' OutFastqFile ResolutionTsv')
    fastq_file = generate_fastq()
    primers = list(SeqIO.parse(fastq_file, 'fastq'))
    consensus_len = len(primers[0])
    name = argv[1].split('.')[0]
    reverse = primers[2].reverse_complement(id=primers[2].id)
    forward = primers[1]
    forward_info = get_id_info(forward)
    forward_location = forward_info[-1]
    forward_seq = str(forward.seq).replace('-', '')
    reverse_info = get_id_info(reverse)
    reverse_location = reverse_info[-1]
    reverse_seq = str(reverse.seq).replace('-', '')
    average_primer_length = (len(forward_seq) + len(reverse_seq)) / 2
    product_len_without_primer = abs(reverse_location -
                                     forward_location) - average_primer_length
    resolution = get_resolution(forward_location, reverse_location)
    coverage = min(forward_info[3], reverse_info[3])
    with open('primer_info.tsv', 'a') as info:
        info.write('Name\tType\tResolution\tAlignmentLength\t'
                   'ProductLength\tPrimerCoverage\t'
                   'Forward\tReverse\n')
        info.write('{}\t{}\t{:.2%}\t{}\t{}\t{:.2%}\t{}\t{}\n'.format(
            name, forward_info[0], resolution, consensus_len,
            product_len_without_primer, coverage, forward_seq, reverse_seq))
        info.write('\n')


if __name__ == '__main__':
    main()
