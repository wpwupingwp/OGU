#!/usr/bin/python3

from Bio import SeqIO
from sys import argv


def generate_fastq():
    seq_file, qual_file, fastq_file = argv[1:]
    with open(seq_file, 'r') as seq, open(qual_file, 'r') as qual:
        fastq = SeqIO.QualityIO.PairedFastaQualIterator(seq, qual)
        SeqIO.write(fastq, fastq_file, 'fastq')
    return fastq_file


def get_id_info(record):
    raw_id = record.id
    print(raw_id)
    name2, tm, cov, sum_bitscore_raw, location = raw_id.split('-')
    return [name2, tm, float(cov[:-2])/100, sum_bitscore_raw, float(location)]


def main():
    fastq_file = generate_fastq()
    primers = list(SeqIO.parse(fastq_file, 'fastq'))
    consensus_len = len(primers[0])
    name = argv[1].split('.')[0]
    if len(primers) == 3:
        reverse = primers[2].reverse_complement(id=primers[2].id)
        have_mid_primer = False
    else:
        reverse = primers[-1].reverse_complement(id=primers[-1].id)
        have_mid_primer = True
    forward = primers[1]
    forward_info = get_id_info(forward)
    forward_seq = str(forward.seq).replace('-', '')
    reverse_info = get_id_info(reverse)
    reverse_seq = str(reverse.seq).replace('-', '')
    product_len_without_primer = abs(reverse_info[-1] - forward_info[-1])
    coverage = min(forward_info[2], reverse_info[2])
    with open('primer_info.tsv', 'a') as info:
        info.write('Name\tType\tLen_Alignment(with_up/downstream400bp)\t'
                   'Product_Len\tcov\tForward\tReverse\t\Have_mid_primer\n')
        info.write('\t'.join([name, forward_info[0], str(consensus_len),
                              str(product_len_without_primer), str(coverage),
                              forward_seq, reverse_seq,
                              str(have_mid_primer)]))
        info.write('\n')


if __name__ == '__main__':
    main()
