#!/usr/bin/python3

import argparse
import os
from functools import wraps
from timeit import default_timer as timer
from datetime import datetime


def print_time(function):
    @wraps(function)
    def wrapper(*args, **kargs):
        start = timer()
        result = function(*args, **kargs)
        end = timer()
        print('Cost {1:3f}s.\n'.format(end-start))
        return result
    return wrapper


@print_time
def function():
    pass


def parse_args():
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=main.__doc__)
    arg.add_argument('-taxon', required=True, help='Taxonomy name')
    arg.add_argument('-continue', action='store_true',
                     help='continue broken download process')
    arg.add_argument('-email', default='wpwupingwp@outlook.com',
                     help='email address used by NCBI Genbank')
    arg.add_argument('-group', default='plants',
                     choices=('animals', 'plants', 'fungi', 'protists',
                              'bacteria', 'archaea', 'viruses'),
                     help='Species kind')
    arg.add_argument('-min_len', default=100, type=int, help='minium length')
    arg.add_argument('-max_len', default=10000, type=int,
                     help='maximum length')
    arg.add_argument('-molecular', choices=('DNA', 'RNA'), default='DNA',
                     help='molecular type')
    arg.add_argument('-organelle', default='plastid',
                     choices=('mitochondrion', 'plastid', 'chloroplast'),
                     help='organelle type')
    arg.add_argument('-out',  help='output directory')
    arg.print_help()
    return arg.parse_args()


#biomol_genomic[PROP]
#biomol_mrna[PROP]
#[filter]
def main():
    """Get data from Genbank.
    """
    arg = parse_args()
    # start here
    if arg.out is None:
        arg.out = datetime.now().isoformat().replace(':', '-')
    os.mkdir(arg.out)
    function()
    # end


if __name__ == '__main__':
    main()
