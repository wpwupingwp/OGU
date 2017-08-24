#!/usr/bin/python3

from collections import Counter
from functools import wraps
from timeit import default_timer as timer
from Bio import AlignIO
import argparse
import numpy as np


def print_time(function):
    @wraps(function)
    def wrapper(*args, **kargs):
        start = timer()
        result = function(*args, **kargs)
        end = timer()
        print('The function {0} Cost {1:3f}s.\n'.format(
            function.__name__, end-start))
        return result
    return wrapper


@print_time
def test_format(file_name):
    with open(file_name, 'r') as _:
        line = _.readline()
        if line.startswith('>'):
            return 'fasta'
        elif line.startswith('#'):
            return 'nexus'
        else:
            raise ValueError('Only support fasta and nexus! Quit now.')


@print_time
def convert(old):
    new = np.array([list(i) for i in old])
    return new


@print_time
def read_alignment(input_file, file_format):
    alignment = AlignIO.read(input_file, file_format)
    return alignment, len(alignment), alignment.get_alignment_length()


@print_time
def remove_gap(alignment, length, width):
    useful = 'ATCG'
    empty = '-N'
    # get alignment head
    new = alignment[:, 0:0]
    for index in range(width):
        column = alignment[:, index:(index+1)]
        string = column[:, 0]
        string = string.upper()
        count = Counter(string)
        for letter in empty:
            if count[letter] == length:
                print('Empty in column {}'.format(index))
                break
            elif (count[letter] + 1) == length:
                print('Only one in column {}'.format(index))
                break
            else:
                new = new + column
                break
        for letter in useful:
            if count[letter] == length:
                print('All same in column {}'.format(index))
                break
    return new, new.get_alignment_length()


def parse_args():
    arg = argparse.ArgumentParser(description=main.__doc__)
    arg.add_argument('input', help='input alignment file')
    arg.add_argument('-o', '--output', default='out',
                     help='output directory')
    arg.print_help()
    return arg.parse_args()


def main():
    """docstring
    """
    arg = parse_args()

    file_format = test_format(arg.input)
    alignment, length, width = read_alignment(arg.input, file_format)
    new_alignment, new_width = remove_gap(alignment, length, width)
    AlignIO.write(new_alignment, arg.output, 'fasta')
    print('Remove {} columns.'.format(width-new_width))


if __name__ == '__main__':
    main()
