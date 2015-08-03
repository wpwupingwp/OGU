import sys
from Bio import SeqIO

def main():
    handle = open(sys.argv[1], 'r')
    data = SeqIO.parse(handle, 'gb')
    target = list()
    for record in data:
        organism = record.annotations['organism'].replace(' ', '_')
        accession = record.annotations['accessions'][0]
        full = str(record.seq)
        target.append([
            organism, 
            accession, 
            'full', 
            full
            ])
        for feature in record.features:
            sequence = list()
            position = list()
            if feature.type != 'gene' or 'gene' not in feature.qualifiers: 
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
            for n, frag in enumurate(position):
                sequence = str(record.seq[frag[0]:frag[1]])
                name = str(feature.qualifiers['gene'][0]).replace(' ', '_')
                if n > 0:
                    name = '-'.join([name, n ])
                target.append([organism, accession, name, sequence])

#Output
    for item in target:
        handle = open(item[2]+'.fasta', 'a')
        handle.write(''.join(['>','|'.join([item[0], item[1], item[2]]),'\n',item[3],'\n']))
        handle.close()
    print('Done.')

if __name__ =='__main__':
    main()
