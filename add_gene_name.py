#!/usr/bin/python3

import argparse
import re
import os
from functools import wraps
from glob import glob
from timeit import default_timer as timer


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
def get_gene(file_list):
    file_gene = dict()
    pattern = re.compile(r'>.*_(.*)$')
    for fasta in file_list:
        with open(fasta, 'r') as raw:
            first_line = raw.readline()
            match = re.match(pattern, first_line)
            gene = match.group(1)
            file_gene[fasta] = gene
    return file_gene


@print_time
def rename(file_gene):
    for fasta, gene in file_gene.items():
        os.rename(fasta, '-'.join([gene, fasta]))


def main():
    """
    The fasta file looks like:
    >AB123456_gene
    """
    start_time = timer()
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('--path', default='./',
                        help='target path, default is "./"')
    parser.print_help()
    arg = parser.parse_args()
    vars(arg)
    # start here
    file_list = glob('*.fa*')
    file_gene = get_gene(file_list)
    rename(file_gene)
    # end
    end_time = timer()
    print('Cost {:.3f}s.\n'.format(end_time-start_time))


if __name__ == '__main__':
    main()
