#!/usr/bin/python3

import argparse
import os
from Bio import SeqIO
from timeit import default_timer as timer


def divide_barcode():
    # get barcode dict
    barcode = dict()
    with open(arg.barcode_file, 'r') as input_file:
        for line in input_file:
            if line.startswith(','):
                continue
            line = line.split(sep=',')
            barcode[line[0]] = line[1].strip()
    # analyze input files
    for record in SeqIO.parse(arg.input, 'fastq'):
        if record.seq[:8] in barcode:
            filename = os.path.join(
                arg.out,
                '{0}-{1}.fastq'.format(os.path.splitext(arg.input)[0],
                                       barcode[record.seq[:8]]))
            with open(filename, 'a') as output:
                SeqIO.write(record[30:-30], output, 'fastq')


def main():
    start_time = timer()
    parser = argparse.ArgumentParser()
    parser.add_argument('--barcode_length', dest='barcode_len', default=10,
                        type=int, help='length of barcode')
    parser.add_argument('-b', dest='barcode_file',
                        help='csv file containing barcode info')
    parser.add_argument('input', help='input file, fastq format')
    parser.add_argument('--cut', default=30, type=int, help='cut head and tail')
    parser.add_argument('-o', dest='out', default='out', help='output path')
    global arg
    arg = parser.parse_args()
    if not os.path.exists(arg.out):
        os.mkdir(arg.out)
    divide_barcode()
    end_time = timer()
    print('Finished with {0:.3f}s. You can find results in {1}.\n'.format(
        end_time-start_time, arg.out))


if __name__ == '__main__':
    main()
