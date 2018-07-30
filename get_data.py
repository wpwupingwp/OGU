#!/usr/bin/python3

import argparse
import json
import re
from datetime import datetime
from os import mkdir, sched_getaffinity
from os.path import join as join_path
from subprocess import run
from timeit import default_timer as timer
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation


def parse_args():
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=main.__doc__)
    arg.add_argument("-query", help='query text')
    arg.add_argument('-continue', action='store_true',
                     help='continue broken download process')
    arg.add_argument('-email', default='',
                     help='email address used by NCBI Genbank')
    output = arg.add_argument_group('output')
    output.add_argument('-out',  help='output directory')
    output.add_argument('-rename', action='store_true',
                        help='try to rename gene')
    output.add_argument('-no_expand', default=False,
                        action='store_true',
                        help='do not expand upstream/downstream')
    output.add_argument('-expand_n', type=int, default=200,
                        help='expand length')
    filters = arg.add_argument_group('filters')
    filters.add_argument('-group',
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
    parsed = arg.parse_args()
    if parsed.organelle is not None:
        # 10k to 1m seems enough
        parsed.min_len = 10000
        parsed.max_len = 1000000
    if parsed.query is None and parsed.taxon is None:
        arg.print_help()
        raise ValueError('Please give more specific query!')
    return parsed


def check_tools():
    print('Checking dependencies ..')
    for tools in ('mafft', 'iqtree'):
        check = run('{} --version'.format(tools), shell=True)
        if check.returncode != 0:
            raise Exception('{} not found! Please install it!'.format(tools))
        print('{} OK.'.format(tools))
    check = run('blastn -version', shell=True)
    if check.returncode != 0:
        raise Exception('BLAST not found! Please install it!')
    print('BLAST OK.\n')


def get_query_string(arg):
    condition = list()
    condition.append('("{}"[SLEN] : "{}"[SLEN])'.format(arg.min_len,
                                                        arg.max_len))
    if arg.group is not None:
        condition.append('{}[filter]'.format(arg.group))
    if arg.query is not None:
        condition.append('"{}"'.format(arg.query))
    if arg.molecular is not None:
        d = {'DNA': 'biomol_genomic[PROP]',
             'RNA': 'biomol_mrna[PROP]'}
        condition.append(d[arg.molecular])
    if arg.taxon is not None:
        condition.append('"{}"[ORGANISM]'.format(arg.taxon))
    if arg.organelle is not None:
        condition.append('{}[filter]'.format(arg.organelle))
    return ' AND '.join(condition)


def download(arg, query):
    print('Your query:')
    print(query)
    Entrez.email = arg.email
    query_handle = Entrez.read(Entrez.esearch(db='nuccore', term=query,
                                              usehistory='y'))
    count = int(query_handle['Count'])
    print('{} records found.'.format(count))
    print('Downloading... Ctrl+C to quit')
    json_file = join_path(arg.out, 'query.json')
    with open(json_file, 'w') as _:
        json.dump(query_handle, _)

    if arg.query is None:
        name = safe(arg.taxon)
    else:
        name = safe(arg.query)
    file_name = join_path(arg.out, name+'.gb')
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


def write_seq(name, sequence_id, feature, whole_seq, path, arg):
    """
    Write fasta file.
    """
    filename = join_path(path, name+'.fasta')
    sequence = feature.extract(whole_seq)

    with open(filename, 'a') as handle:
        handle.write(sequence_id+'\n')
        handle.write(str(sequence)+'\n')
    if not arg.no_expand:
        if feature.location_operator == 'join':
            loc = feature.location.parts
            # ensure increasing order
            # parts do not have sort method
            sorted(loc, key=lambda x: x.start)
            new_loc = sum([
                # avoid IndexError
                FeatureLocation(max(0, loc[0].start-arg.expand_n),
                                loc[0].end, loc[0].strand),
                *loc[1:-1],
                FeatureLocation(loc[-1].start,
                                min(len(whole_seq), loc[-1].end+arg.expand_n),
                                loc[-1].strand)])
            feature.location = new_loc
        feature.type = 'expand'
        sequence = feature.extract(whole_seq)
        filename2 = join_path(path, 'expand.{}.fasta'.format(name))
        with open(filename2, 'a') as handle:
            handle.write(sequence_id+'\n')
            handle.write(str(sequence)+'\n')
        return filename2
    return filename


def get_feature_name(feature, arg):
    """
    Get feature name and collect genes for extract spacer.
    Only handle gene, product, misc_feature, misc_RNA.
    Return: [name, feature.type]
    """
    name = None
    if feature.type == 'gene':
        if 'gene' in feature.qualifiers:
            gene = feature.qualifiers['gene'][0].replace(' ', '_')
            if arg.rename:
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
        if (('intergenic_spacer' in misc_feature or
             'IGS' in misc_feature)):
            # 'IGS' in misc_feature) and len(misc_feature) < 100):
            name = safe(misc_feature)
            name = name.replace('intergenic_spacer_region',
                                'intergenic_spacer')
        else:
            pass
    elif feature.type == 'misc_RNA':
        if 'product' in feature.qualifiers:
            misc_feature = feature.qualifiers['product'][0].replace(
                ' ', '_')
        elif 'note' in feature.qualifiers:
            misc_feature = feature.qualifiers['note'][0].replace(
                ' ', '_')
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
        pass
        # print('Cannot handle feature:')
        # print(feature)
    return name, feature.type


def get_spacer(genes, arg):
    """
    List: [[name, SeqFeature],]
    """
    spacers = list()
    # sorted according to sequence starting postion
    genes.sort(key=lambda x: int(x[1].location.start))
    for n, present in enumerate(genes[1:], 1):
        before = genes[n-1]
        # use sort to handle complex location relationship of two fragments
        location = [before[1].location.start, before[1].location.end,
                    present[1].location.start, present[1].location.end]
        location.sort(key=lambda x: int(x))
        start, end = location[1:3]
        if before[1].location.strand == present[1].location.strand == -1:
            strand = -1
        else:
            strand = 1
        name = '_'.join([before[0], present[0]])
        spacer = SeqFeature(FeatureLocation(start, end), id=name,
                            type='spacer', strand=strand)
        spacers.append(spacer)
    return spacers


def divide(gbfile, arg):
    """
    Given genbank file, return divided fasta files.
    """
    start = timer()
    groupby_gene = '{}-groupby_gene'.format(gbfile.replace('.gb', ''))
    mkdir(groupby_gene)
    groupby_name = '{}-groupby_name'.format(gbfile.replace('.gb', ''))
    mkdir(groupby_name)
    handle_raw = open(gbfile+'.fasta', 'w')
    wrote_by_gene = set()
    wrote_by_name = set()

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
        genes = list()

        for feature in record.features:
            name, feature_type = get_feature_name(feature, arg)
            # skip unsupport feature
            if name is None:
                continue
            # skip abnormal annotation
            if len(feature) > 20000:
                print('Skip abnormal annotaion of {}!'.format(name))
                print('Accession: ', accession)
                continue
            if feature_type == 'gene':
                genes.append([name, feature])
            feature_name.append(name)
            sequence_id = '>' + '|'.join([name, taxon, accession, specimen])
            wrote = write_seq(name, sequence_id, feature, whole_seq,
                              groupby_gene, arg)
            wrote_by_gene.add(wrote)

        # extract spacer
        spacers = get_spacer(genes, arg)
        for spacer in spacers:
            sequence_id = '>' + '|'.join([spacer.id, taxon,
                                          accession, specimen])
            wrote = write_seq(spacer.id, sequence_id, spacer, whole_seq,
                              groupby_gene, arg)
            wrote_by_gene.add(wrote)
        # write to group_by name, i.e., one gb record one fasta
        if 'ITS' in feature_name:
            name_str = 'ITS'
        elif len(feature_name) >= 4:
            name_str = '{}-...-{}'.format(feature_name[0], feature_name[-1])
        elif len(feature_name) == 0:
            name_str = 'Unknown'
        else:
            name_str = '-'.join(feature_name)
        # directly use genome type as name
        if arg.organelle is not None:
            name_str = '{}_genome'.format(arg.organelle)
        record.id = '|'.join([name_str, taxon, accession, specimen])
        record.description = ''
        filename = join_path(groupby_name, name_str+'.fasta')
        with open(filename, 'a') as out:
            SeqIO.write(record, out, 'fasta')
            wrote_by_name.add(filename)
        # write raw fasta
        SeqIO.write(record, handle_raw, 'fasta')

    end = timer()
    print('Divide done with {:.3f}s.'.format(end-start))
    return wrote_by_gene, wrote_by_name


def mafft(files):
    result = list()
    # get available CPU cores
    cores = len(sched_getaffinity(0))
    print('Start mafft ...')
    for fasta in files:
        print('Aligning {}'.format(fasta))
        out = fasta + '.aln'
        _ = ('mafft --thread {} --reorder --quiet --adjustdirection '
             ' {} > {}'.format(cores-1, fasta, out))
        m = run(_, shell=True)
        if m.returncode == 0:
            result.append(out)
    print('Done with mafft.')
    return result


def main():
    """Get data from Genbank.
    """
    arg = parse_args()
    if arg.out is None:
        arg.out = datetime.now().isoformat().replace(':', '-')
    mkdir(arg.out)
    check_tools()
    query = get_query_string(arg)
    gbfile = download(arg, query)
    wrote_by_gene, wrote_by_name = divide(gbfile, arg)
    if arg.max_len > 10000:
        # two few sequences will cause empty output, which was omit
        aligned = mafft(wrote_by_gene)
    else:
        aligned = mafft(wrote_by_name)
    for aln in aligned:
        print(aln)


if __name__ == '__main__':
    main()
