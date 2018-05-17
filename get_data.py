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
        print('Cost {0:3f}s.\n'.format(end-start))
        return result
    return wrapper


@print_time
def function():
    pass


def parse_args():
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=main.__doc__)
    arg.add_argument('query', help='query text')
    arg.add_argument('-continue', action='store_true',
                     help='continue broken download process')
    arg.add_argument('-email', default='wpwupingwp@outlook.com',
                     help='email address used by NCBI Genbank')
    arg.add_argument('-out',  help='output directory')
    filters = arg.add_argument_group('filters')
    filters.add_argument('-group', default='plants',
                         choices=('animals', 'plants', 'fungi', 'protists',
                                  'bacteria', 'archaea', 'viruses'),
                         help='Species kind')
    filters.add_argument('-min_len', default=100, type=int,
                         help='minium length')
    filters.add_argument('-max_len', default=10000, type=int,
                         help='maximum length')
    filters.add_argument('-molecular', choices=('DNA', 'RNA'),
                         help='molecular type')
    filters.add_argument('-taxon', help='Taxonomy name')
    filters.add_argument('-organelle',
                         choices=('mitochondrion', 'plastid', 'chloroplast'),
                         help='organelle type')
    arg.print_help()
    return arg.parse_args()


def get_query_string(arg):
    condition = list()
    condition.append('"{}"'.format(arg.query))
    condition.append('{}[filter]'.format(arg.group))
    condition.append('("{}"[SLEN] : "{}"[SLEN])'.format(arg.min_len,
                                                        arg.max_len))
    if arg.molecular is not None:
        d = {'DNA': 'biomol_genomic[PROP]',
             'RNA': 'biomol_mrna[PROP]'}
        condition.append(d[arg.molecular])
    if arg.taxon is not None:
        condition.append('"{}"[ORGANISM]'.format(arg.taxon))
    if arg.organelle is not None:
        condition.append('{}[filter]'.format(arg.organelle))
    return ' AND '.join(condition)


def main():
    """Get data from Genbank.
    """
    arg = parse_args()
    # start here
    if arg.out is None:
        arg.out = datetime.now().isoformat().replace(':', '-')
    query = get_query_string(arg)
    print(query)
    os.mkdir(arg.out)
    function()
    # end


if __name__ == '__main__':
    main()
