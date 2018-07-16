#!/usr/bin/python3

import argparse
import json
import re
from datetime import datetime
from os import mkdir
from os.path import join as join_path
from subprocess import run
from timeit import default_timer as timer
from Bio import Entrez, SeqIO
from Bio.Seq import Seq


def check_tools():
    for tools in ('mafft', 'iqtree'):
        check = run('{} --veriosn'.format(tools), shell=True)
        if check.returncode != 0:
            raise Exception('{} not found! Please install it!'.format(tools))
        check = run('blastn -version', shell=True)
        if check.returncode != 0:
            raise Exception('BLAST not found! Please install it!')


def download(arg, query):
    Entrez.email = arg.email
    query_handle = Entrez.read(Entrez.esearch(db='nuccore', term=query,
                                              usehistory='y'))
    count = int(query_handle['Count'])
    print('Your query:')
    print(query)
    print('{} records found.'.format(count))
    print('Downloading... Ctrl+C to quit')
    json_file = join_path(arg.out, 'query.json')
    with open(json_file, 'w') as _:
        json.dump(query_handle, _)

    file_name = join_path(arg.out, arg.query+'.gb')
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


def gene_rename(old_name):
    """For chloroplast gene.
    Input->str
    Output->List[new_name:str, name_type:str]
    """
    lower = old_name.lower()
    # (trna|trn(?=[b-z]))
    s = re.compile(r'(\d+\.?\d?)(s|rrn|rdna)')
    if lower.startswith('trn'):
        pattern = re.compile(r'([atcgu]{3})')
        try:
            codon = Seq(re.search(pattern, lower).group(1))
        except AttributeError:
            return old_name, 'bad_name'
        try:
            new_name = 'trn{}{}'.format(codon.reverse_complement().translate(),
                                        codon.transcribe())
        except ValueError:
            return old_name, 'bad_name'
        gene_type = 'tRNA'
    elif lower.startswith('rrn'):
        pattern = re.compile(r'(\d+\.?\d?)')
        try:
            number = re.search(pattern, lower).group(1)
        except AttributeError:
            return old_name, 'bad_name'
        new_name = 'rrn{}'.format(number)
        gene_type = 'rRNA'
    elif re.search(s, lower) is not None:
        new_name = 'rrn{}'.format(re.search(s, lower).group(1))
        gene_type = 'rRNA'
    else:
        pattern = re.compile(r'[^a-z]*'
                             '(?P<gene>[a-z]+)'
                             '[^a-z0-9]*'
                             '(?P<suffix>[a-z]|[0-9]+)')
        match = re.search(pattern, lower)
        try:
            gene = match.group('gene')
            suffix = match.group('suffix')
        except ValueError:
            return old_name, 'bad_name'
        new_name = '{}{}'.format(gene, suffix.upper())
        # captitalize last letter
        if len(new_name) > 3:
            s = list(new_name)
            if s[-1].isalpha():
                new_name = '{}{}'.format(
                    ''.join(s[:-1]), ''.join(s[-1]).upper())
        gene_type = 'normal'
    if len(lower) >= 15:
        gene_type = 'suspicious_name'
    return new_name, gene_type


def safe(old):
        return re.sub(r'\W', '_', old)


def get_taxon(order_family):
    """
    From Zhang guojin
    order end with ales
    family end with aceae except 8
    http://duocet.ibiodiversity.net/index.php?title=%E4%BA%92%E7%94%A8%E5%90%8D
    %E7%A7%B0&mobileaction=toggle_view_mobile"""
    # order|family|organims(genus|species)
    family_exception_raw = (
        'Umbelliferae,Palmae,Compositae,Cruciferae,Guttiferae,Leguminosae,'
        'Leguminosae,Papilionaceae,Labiatae,Gramineae')
    family_exception = family_exception_raw[0].split(',')
    order = ''
    family = ''
    for item in order_family:
        if item.endswith('ales'):
            order = item
        elif item.endswith('aceae'):
            family = item
        elif item in family_exception:
            family = item
    return order, family


def get_seq(feature, whole_sequence, expand=False, expand_n=100):
    """
    Given location and whole sequence, return fragment.
    """
    if expand:
        # in case of negative start
        start = max(0, feature.location.start-expand_n)
        end = min(len(whole_sequence), feature.location.end+expand_n)
    else:
        start = feature.location.start
        end = feature.location.end
    return whole_sequence[start:end]


def divide(gbfile, rename=True, expand=True):
    start = timer()

    groupby_gene = '{}-groupby_gene'.format(gbfile.replace('.gb', ''))
    mkdir(groupby_gene)
    groupby_name = '{}-groupby_name'.format(gbfile.replace('.gb', ''))
    mkdir(groupby_name)
    handle_raw = open(gbfile+'.fasta', 'w')

    def extract(feature, name, whole_seq):
        filename = join_path(groupby_gene, name+'.fasta')
        sequence = get_seq(feature, whole_seq, expand=False)
        with open(filename, 'a') as handle:
            handle.write('>{}|{}|{}|{}\n{}\n'.format(
                name, taxon, accession, specimen, sequence))
        filename2 = join_path(groupby_gene, 'expand.{}.fasta'.format(name))
        sequence = get_seq(feature, whole_seq, expand=True, expand_n=100)
        with open(filename2, 'a') as handle:
            handle.write('>{}|{}|{}|{}\n{}\n'.format(
                name, taxon, accession, specimen, sequence))

    for record in SeqIO.parse(gbfile, 'gb'):
        # only accept gene, product, and spacer in misc_features.note
        order_family = record.annotations['taxonomy']
        order, family = get_taxon(order_family)
        organism = record.annotations['organism'].replace(' ', '_')
        genus, *species = organism.split('_')
        taxon = '{}|{}|{}|{}'.format(order, family, genus, '_'.join(species))
        accession = record.annotations['accessions'][0]
        try:
            specimen = record.features[0].qualifiers['specimen_voucher'
                                                     ][0].replace(' ', '_')
        except (IndexError, KeyError):
            specimen = ''
        whole_seq = record.seq
        feature_name = list()

        for feature in record.features:
            gene = ''
            if feature.type == 'gene':
                if 'gene' in feature.qualifiers:
                    gene = feature.qualifiers['gene'][0].replace(' ', '_')
                    gene = gene_rename(gene)[0]
                    name = safe(gene)
                elif 'product' in feature.qualifiers:
                    product = feature.qualifiers['product'][0].replace(
                        ' ', '_')
                    name = safe(product)
            elif feature.type == 'misc_feature':
                if 'product' in feature.qualifiers:
                    misc_feature = feature.qualifiers['product'][0].replace(
                        ' ', '_')
                elif 'note' in feature.qualifiers:
                    misc_feature = feature.qualifiers['note'][0].replace(
                        ' ', '_')
                else:
                    continue
                if (('intergenic_spacer' in misc_feature or
                     'IGS' in misc_feature) and len(misc_feature) < 100):
                    name = safe(misc_feature)
                    name = name.replace('intergenic_spacer_region',
                                        'intergenic_spacer')
                else:
                    print('Too long name: {}'.format(misc_feature))
            elif feature.type == 'misc_RNA':
                if 'product' in feature.qualifiers:
                    misc_feature = feature.qualifiers['product'][0].replace(
                        ' ', '_')
                elif 'note' in feature.qualifiers:
                    misc_feature = feature.qualifiers['note'][0].replace(
                        ' ', '_')
                else:
                    continue
                name = safe(misc_feature)
                # handle ITS
                if 'internal_transcribed_spacer' in name:
                    name = 'ITS'
                # name = name.replace('internal_transcribed_spacer', 'ITS')
                # if 'ITS_1' in name:
                #     if 'ITS_2' in name:
                #         name = 'ITS'
                #     else:
                #         name = 'ITS_1'
                # elif 'ITS_2' in name:
                #     name = 'ITS_2'
            else:
                print(feature)
                continue
            extract(feature, name, whole_seq)
            feature_name.append(name)
        if 'ITS' in feature_name:
            name_str = 'ITS'
        elif len(feature_name) >= 4:
            name_str = '{}-{}features-{}'.format(
                feature_name[0], len(feature_name)-2, feature_name[-1])
        elif len(feature_name) == 0:
            name_str = 'Unknown'
        else:
            name_str = '-'.join(feature_name)

        record.id = '|'.join([name_str, taxon, accession, specimen])
        record.description = ''
        SeqIO.write(record, handle_raw, 'fasta')
        with open(join_path(groupby_name, name_str+'.fasta'),
                  'a') as handle_name:
            SeqIO.write(record, handle_name, 'fasta')

    end = timer()
    print('Divide done with {:.3f}s.'.format(end-start))
    return groupby_gene, groupby_name


def parse_args():
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=main.__doc__)
    arg.add_argument('query', help='query text')
    arg.add_argument('-continue', action='store_true',
                     help='continue broken download process')
    arg.add_argument('-email', default='',
                     help='email address used by NCBI Genbank')
    output = arg.add_argument_group('output')
    output.add_argument('-out',  help='output directory')
    output.add_argument('-rename', action='store_true',
                        help='try to rename gene')
    output.add_argument('--expand', type=int, default=200,
                        help='expand length of upstream/downstream')
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
    check_tools()
    arg = parse_args()
    if arg.out is None:
        arg.out = datetime.now().isoformat().replace(':', '-')
    query = get_query_string(arg)
    mkdir(arg.out)
    gbfile = download(arg, query)
    groupby_gene, groupby_name = divide(gbfile, arg.rename, arg.expand)


if __name__ == '__main__':
    main()
