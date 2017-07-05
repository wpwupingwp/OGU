#!/usr/bin/python3

import argparse
import os
from timeit import default_timer as timer
from Bio import SeqIO


def get_format(filename):
    with open(filename, 'r') as raw:
        line = raw.readline()
        if line.startswith('>'):
            return 'fasta'
        elif line.startswith('@'):
            return 'fastq'
        else:
            raise ValueError('Unsupport format!')


def split(fasta, size, out, file_format):
    raw = SeqIO.parse(fasta, file_format)

    def newfile(index):
        output = open(os.path.join(
            out, '{}.{}'.format(fasta, (index//size)+1)), 'w')
        return output

    for index, record in enumerate(raw, 1):
        if index % size == 1:
            output = newfile(index)
        SeqIO.write(record, output, file_format)
    return


def parse_args():
    arg = argparse.ArgumentParser(description=main.__doc__)
    arg.add_argument('-i', '--input', required=True, help='input file')
    arg.add_argument('-s', '--size', type=int, default=100000,
                     help='how many sequences one file have')
    arg.add_argument('-o', '--out', help='output directory')
    arg.print_help()
    return arg.parse_args()


def main():
    """docstring
    """
    start = timer()
    arg = parse_args()

    file_format = get_format(arg.input)
    if arg.out is None:
        arg.out = 'out_{}'.format(arg.input)
    os.mkdir(arg.out)
    split(arg.input, arg.size, arg.out, file_format)
    # end
    end = timer()
    print('Cost {:3f}s.\n'.format(end-start))


if __name__ == '__main__':
    main()
