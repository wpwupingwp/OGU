#!/usr/bin/python3

import sys
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord 
from Bio.Blast.Applications import NcbiblastnCommandline as nb
from os import makedirs
from os.path import exists

def get_cds(fragments, option):
    if option == '1':
        pass
    elif option == '2':
        return 
    else:
        raise ValueError('Wrong option!\n')

    wanted_gene = [
        'accD', 'atpA', 'atpB', 'atpE', 'atpF', 'atpH', 'atpI', 'ccsA', 'cemA', 
        'clpP', 'infA', 'matK', 'ndhA', 'ndhB', 'ndhC', 'ndhD', 'ndhE', 'ndhF', 
        'ndhG', 'ndhH', 'ndhI', 'ndhJ', 'ndhK', 'petA', 'petB', 'petD', 'petG', 
        'petL', 'petN', 'psaA', 'psaB', 'psaC', 'psaI', 'psaJ', 'psbA', 'psbB', 
        'psbC', 'psbD', 'psbE', 'psbF', 'psbH', 'psbI', 'psbJ', 'psbK', 'psbL', 
        'psbM', 'psbN', 'psbT', 'psbZ', 'rbcL', 'rpl14', 'rpl16', 'rpl2', 
        'rpl20', 'rpl22', 'rpl23', 'rpl32', 'rpl33', 'rpl36', 'rpoA', 'rpoB', 
        'rpoC1', 'rpoC2', 'rps11', 'rps12', 'rps14', 'rps15', 'rps16', 'rps18', 
        'rps19', 'rps2', 'rps3', 'rps4', 'rps7', 'rps8', 'rrn16', 'rrn23', 
        'rrn4.5', 'rrn5', 'ycf1', 'ycf2', 'ycf3', 'ycf4'
    ]
    handle = open(sys.argv[1], 'r')
    data = SeqIO.parse(handle, 'gb')
    for record in data:
        organism = record.annotations['organism'].replace(' ', '_')
        accession = record.annotations['accessions'][0]
        for feature in record.features:
            sequence = list()
            position = list()
            if feature.type != 'CDS' or 'gene' not in feature.qualifiers: 
                continue
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
                fragments.append([organism, accession, name, sequence])

def out_cds(fragments, option):
    if option == '1':
        pass
    else:
        return
    handle_all = open('output/all.fasta', 'a')
    for item in fragments:
        handle = open(''.join([
            'output/',
            item[2],
            '.fasta'
        ]),
                      'a')
        handle.write(''.join(['>','|'.join([item[0], item[1], item[2]]),'\n',item[3],'\n']))
        handle_all.write(''.join(['>','|'.join([item[0], item[1], item[2]]),'\n',item[3],'\n']))
        handle.close()
    handle_all.close()
    print('Done.')

def blast(option):
    if option == '1':
        query_file = './output/all.fasta'
    else:
        query_file = sys.argv[1]
    cmd = nb(
    #  num_threads=8,
        query=query_file,
        db=sys.argv[2], 
        task='megablast', 
        evalue=0.001, 
        outfmt=5, 
        # xml format
        out='BlastResult.xml'
    )
    stdout, stderr = cmd()
    return 

def parse(target):
    blast_result = list(SearchIO.parse('BlastResult.xml', 'blast-xml'))
    for record in blast_result:
        if len(record) == 0:
            continue
        else:
            tophit = record[0]
        target.append([tophit[0][0].query, tophit[0][0].hit])

def output(target, option):
    if option == '1':
        for record in target:
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
    else:
        output_file = open('output/' + sys.argv[1] + '-filtered.fasta', 'w' )
        contig_id = {i[0].id for i in target} 
        query_file = SeqIO.parse(sys.argv[1], 'fasta')
        for record in query_file:
            if record.id in contig_id:
                SeqIO.write(record, output_file, 'fasta')

def main():
    """Usage:
    python3 getCDS.py genbank_file/fasta_file contig_file mode
    Before running this script, ensure you already put blast database file in
    current path. Also, it assumes that you have installed BLAST suite before. 
    To create the db file: 
    makeblastdb -in infile -out outfile -dbtype nucl
    Before run, you need to chose mode:
        1. Query CDS in contig
        2. Query contig in genome"""
    if exists('output') == False:
        makedirs('output')
    fragments = list()
    target = list()
    print('\n', main.__doc__)
    option = sys.argv[3]
    get_cds(fragments, option)
    out_cds(fragments, option)
    blast(option)
    parse(target)
    output(target, option)

if __name__ =='__main__':
    main()
