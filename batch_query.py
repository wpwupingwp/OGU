#!/usr/bin/python3

from Bio import Entrez
import argparse
import os


def get_list(list_file):
    down_list = list()
    with open(list_file, 'r') as raw:
        for line in raw:
            down_list.append(line.strip())
    return down_list


def down(taxon_name):
    EMAIL = 'wpwupingwp@outlook.com'
    Entrez.email = EMAIL
    FILTER = 'plastid'
    # set 'noexp' to fetch only this level
    EXP = 'exp'
    query = '''{}[Organism:{}] AND ({}[filter] AND ("{}"[SLEN] :
    "{}"[SLEN]))'''.format(taxon_name, EXP, FILTER, arg.min_len, arg.max_len)
    print(query)
    handle = Entrez.read(Entrez.esearch(db='nuccore', term=query,
                                        usehistory='y',))
    genome_content = Entrez.efetch(db='nuccore',
                                   webenv=handle['WebEnv'],
                                   query_key=handle['QueryKey'],
                                   rettype='gb',
                                   retmode='text')
    file_name = os.path.join(arg.out, taxon_name+'.gb')
    with open(file_name, 'w') as output_file:
            output_file.write(genome_content.read())


def main():
    global arg
    arg = argparse.ArgumentParser()
    arg.add_argument('-l', dest='list_file',
                     help='Taxonomy name list, seperate by line')
    arg.add_argument('-o', dest='out', default='out', help='output path')
    arg.add_argument('-min_len', default=10, type=int, help='minium length')
    arg.add_argument('-max_len', default=200000, type=int,
                     help='maximum length')
    arg = arg.parse_args()
    down_list = get_list(arg.list_file)
    down(down_list[1])


if __name__ == '__main__':
    main()
