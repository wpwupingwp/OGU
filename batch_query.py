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

def down(down_list):
    EMAIL = 'wpwupingwp@outlook.com'
    Entrez.email = EMAIL
    query = '{0}[Organism:exp] AND (plastid[filter] AND ("10000"[SLEN] : "50000"[SLEN]))'.format('txid4479')
    handle = Entrez.read(Entrez.esearch(db='nuccore', term=query,
                                        usehistory='y',))
    genome_content = Entrez.efetch(db='nuccore',
                                   webenv=handle['WebEnv'],
                                   query_key=handle['QueryKey'],
                                   rettype='gb',
                                   retmode='text')
    with open('test', 'w') as output_file:
            output_file.write(genome_content.read())


def main():
    arg = argparse.ArgumentParser()
    arg.add_argument('list_file', help='Taxonomy name list, seperate by line')
    arg.add_argument('-o', dest='output', default='out', help='output path')
    arg.add_argument('-min_len', default=10, type=int, help='minium length')
    arg.add_argument('-max_len', default=50000, type=int,
                     help='maximum length')
