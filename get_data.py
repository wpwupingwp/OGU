#!/usr/bin/python3

import argparse
import json
import os
from datetime import datetime
from Bio import Entrez


def download(arg, query):
    Entrez.email = arg.email
    query_handle = Entrez.read(Entrez.esearch(db='nuccore', term=query,
                                              usehistory='y'))
    count = int(query_handle['Count'])
    print('Your query:')
    print(query)
    print('{} records found.'.format(count))
    print('Downloading... Ctrl+C to quit')
    json_file = os.path.join(arg.out, 'query.json')
    with open(json_file, 'w') as _:
        json.dump(query_handle, _)

    file_name = os.path.join(arg.out, arg.query+'.gb')
    output = open(file_name, 'w')
    ret_start = 0
    ret_max = 1000
    while ret_start <= count:
        print('{}-{}'.format(ret_start, ret_start+ret_max))
        try:
            data = Entrez.efetch(db='nuccore',
                                 webenv=query_handle['WebEnv'],
                                 query_key=query_handle['QueryKey'],
                                 rettype='gb',
                                 retmode='text',
                                 retstart=ret_start,
                                 retmax=ret_max)
            output.write(data.read())
        # just retry if connection failed
        except IOError:
            print('Retrying...')
            continue
        ret_start += 1000
    print('Done.')
    return file_name


def parse_args():
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=main.__doc__)
    arg.add_argument('query', help='query text')
    arg.add_argument('-continue', action='store_true',
                     help='continue broken download process')
    arg.add_argument('-email', default='',
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
    if arg.out is None:
        arg.out = datetime.now().isoformat().replace(':', '-')
    query = get_query_string(arg)
    os.mkdir(arg.out)
    download(arg, query)


if __name__ == '__main__':
    main()
