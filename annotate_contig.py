#!/usr/bin/python3

import sys
from Bio import SeqIO
from Bio import SearchIO
from Bio.Blast.Applications import NcbiblastnCommandline as nb
from Bio.SeqRecord import SeqRecord 
from multiprocessing import cpu_count
from os import makedirs
from os.path import exists
from subprocess import call


def get_cds():
    #Edit it when necessary
    wanted_gene_list = [
        'accD', 'atpA', 'atpB', 'atpE', 'atpF', 'atpH', 'atpI', 'ccsA',
        'cemA', 'clpP', 'infA', 'matK', 'ndhA', 'ndhB', 'ndhC', 'ndhD',
        'ndhE', 'ndhF', 'ndhG', 'ndhH', 'ndhI', 'ndhJ', 'ndhK', 'petA',
        'petB', 'petD', 'petG', 'petL', 'petN', 'psaA', 'psaB', 'psaC',
        'psaI', 'psaJ', 'psbA', 'psbB', 'psbC', 'psbD', 'psbE', 'psbF',
        'psbH', 'psbI', 'psbJ', 'psbK', 'psbL', 'psbM', 'psbN', 'psbT',
        'psbZ', 'rbcL', 'rpl14', 'rpl16', 'rpl2', 'rpl20', 'rpl22',
        'rpl23', 'rpl32', 'rpl33', 'rpl36', 'rpoA', 'rpoB', 'rpoC1',
        'rpoC2', 'rps11', 'rps12', 'rps14', 'rps15', 'rps16', 'rps18',
        'rps19', 'rps2', 'rps3', 'rps4', 'rps7', 'rps8', 'rrn16',
        'rrn23', 'rrn4.5', 'rrn5', 'ycf1', 'ycf2', 'ycf3', 'ycf4'
    ]
    cds = list()
    genomes = SeqIO.parse(sys.argv[1], 'gb')
    for genome in genomes:
        for feature in genome.features:
            if feature.type != 'CDS' or 'gene' not in feature.qualifiers: 
                continue
            sequence = list()
            position = list()
            if feature.location_operator != 'join':
                position.append([
                    int(feature.location.start), 
                    int(feature.location.end)
                ])
            else:
                for i in feature.sub_features:
                    position.append([
                        int(i.location.start), 
                        int(i.location.end)
                    ])
            for n, frag in enumerate(position):
                name = str(feature.qualifiers['gene'][0]).replace(' ', '_')
                if name not in wanted_gene:
                    continue
                sequence = str(record.seq[frag[0]:frag[1]])
                if n > 0:
                    name = '-'.join([name, str(n+1)])
                cds.append([name, sequence])
    return cds

def out_cds(cds):
    handle = open('cds.fasta' ,'w')
    for gene in cds:
        handle.write(''.join([
            '>',gene[0], '|', gene[1], '\n',
            gene[2],'\n'
        ]))
    handle.close()
    return 'cds.fasta'


def blast(query_file):
    contig_file = sys.argv[2]
    xml_file = 'out/BlastResult.xml'
    call('makeblastdb -in {0} -out {1} -dbtype nucl'.format(contig_file,
                                                            contig_file),
         shell=True)
    cmd = nb(
        num_threads=cpu_count(),
        query='query_file',
        db=contig_file,
        task='blastn', 
        evalue=0.001, 
        outfmt=5, 
        out=xml_file
    )
    stdout, stderr = cmd()
    return xml_file


def parse(xml_file):
    parse_result = list()
    blast_result = SearchIO.parse(xml_file, 'blast-xml')
    for record in blast_result:
        if len(record) == 0:
            continue
        else:
            tophit = record[0]
        parse_result.append([tophit[0][0].query, tophit[0][0].hit])
    return parse_result

def output(result):
    for record in result:
        gene = record[0].id.split(sep='|')[-1]
        output_file = ''.join([
            'output/',
            sys.argv[2], 
            '-',
            gene, 
            '.fasta'
        ])
        rename_seq = SeqRecord(
            seq=record[1].seq, 
            id='|'.join([
                gene,
                sys.argv[1],
                record[1].id
            ]),
            description=''
        )
        SeqIO.write(rename_seq, output_file, 'fasta')
    :
    output_file = open('output/' + sys.argv[1] + '-filtered.fasta', 'w' )
    contig_id = {i[0].id for i in target} 
    query_file = SeqIO.parse(sys.argv[1], 'fasta')
    for record in query_file:
        if record.id in contig_id:
            SeqIO.write(record, output_file, 'fasta')

def main():
    """
    This program will annotate contigs from assembly according to given
    genbank file, which describes a complete chloroplast genome. The genbank
    file may contains single or several genomes.
    Edit wanted_gene list in get_cds(). If you want to annotate 
    mitochrondria contigs.
    Usage:
    python3 annotate_contig.py genbank_file contig_file mode
    Mode:
        1. Query contig against coding genes, then every contig will be
        annotated by gene name.
        2. Query contig in a whole genome. It only judge if contig was belong to
        genome of given genbank file.
    On default, it use Nicotiana.gb which was placed in current path."""
    print(main.__doc__)
    if not exists('out'):
        makedirs('out')
    if mode not in ['1', '2']:
        raise ValueError('Bad command!\n')
    if mode == '1':
        cds = get_cds()
        query_file = out_cds(cds)
        xml_file = blast(query_file)
        result = parse(xml_file)
        output(result)

if __name__ =='__main__':
    main()
