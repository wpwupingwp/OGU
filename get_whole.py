#!/usr/bin/python3

import argparse
from Bio import SearchIO, SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline as nb
from multiprocessing import cpu_count


def blast(query_file, contig_file):
    """
    BLAST result located in 'out/'."""
    xml_file = 'BlastResult.xml'
    cmd = nb(
        num_threads=cpu_count(),
        query=query_file,
        db=contig_file,
        #  max_hsps=1,
        max_target_seqs=1,
        task='blastn',
        evalue=0.001,
        outfmt=5,
        out=xml_file
    )
    stdout, stderr = cmd()
    return xml_file


def parse(xml_file):
    parse_result = dict()
    blast_result = SearchIO.parse(xml_file, 'blast-xml')
    for record in blast_result:
        if len(record) == 0:
            continue
        for i in record:
            parse_result[i[0][0].hit_id] = i[0][0].query_id
    return parse_result


def output(raw_file, output_file, parse_result):
    handle = open(output_file, 'a')
    raw = SeqIO.parse(raw_file, 'fasta')
    for i in raw:
        new_id = ''
        try:
            new_id = parse_result[i.id]
        except:
            continue
        handle.write('>{0}|{1}\n{2}\n'.format(
            new_id, i.id, i.seq))
    handle.close()


def main():
    parameters = argparse.ArgumentParser(
        description='get whole sequence accoring to BLAST result')
    parameters.add_argument('query', help='fasta file you want to query')
    parameters.add_argument('database', help='BLAST database you use')
    parameters.add_argument('raw', help='''raw fasta file containing all
                            sequence, also used to build BLAST library''')
    parameters.add_argument('-o', '--output',
                            default='output.fasta', help='output file')
    parameters.print_help()
    arg = parameters.parse_args()
    blast_result = blast(arg.query, arg.database)
    parse_result = parse(blast_result)
    output(arg.raw, arg.output, parse_result)


if __name__ == '__main__':
    main()
