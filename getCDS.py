#!/usr/bin/python3
import sys
from Bio import SeqIO

def main():
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
    target = list()
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
                target.append([organism, accession, name, sequence])

#Output
    for item in target:
        handle = open(''.join([
            'output/',
            item[2],
            '.fasta'
        ]),
                      'a')
        handle.write(''.join(['>','|'.join([item[0], item[1], item[2]]),'\n',item[3],'\n']))
        handle.close()
    print('Done.')

if __name__ =='__main__':
    main()
