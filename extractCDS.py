import sys
from Bio import SeqIO

def main():
    handle = open(sys.argv[1], 'r')
    data = SeqIO.parse(handle, 'gb')
    target = list()
    cds = list()
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
#Use for to get all sub_features
        for feature in record.features:
           # try:
          #      position  =  [[i.location.start, i.location.end] for i in feature.sub_features]
          #  except:
          #      position  =  [feature.location.start, feature.location.end]
          #if feature.type == 'misc_feature': 
          #    name = str(feature.qualifiers['note'][0]).replace(' ', '_')
          #elif feature.type =='gene' and 'gene' in feature.qualifiers:
          #    name = str(feature.qualifiers['gene'][0]).replace(' ', '_')

          #if feature.type != 'CDS' and 'CDS' in feature.qualifiers:

            sequence = list()
            position = list()
            if feature.type != 'CDS':
                continue

#some CDS doesn't have name???
            if 'gene' in feature.qualifiers: 
#avoid space by replacing it with dash
                name = str(feature.qualifiers['gene'][0]).replace(' ', '_')
            else:
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
            for frag in position:
                sequence.append(str(record.seq[frag[0]:frag[1]]))

            sequence = ''.join(sequence)
            target.append([organism, accession, name, sequence])

#Output
    for item in target:
        handle = open(item[2], 'a')
        handle.write(''.join(['>','|'.join([item[0], item[1], item[2]]),'\n',item[3],'\n']))
        handle.close()
    #for item in target:
    #    handle_out.write('>')
    #    handle_out.write('|'.join([item[0], item[1],item[2]]))
    #    handle_out.write('\n')
    #    handle_out.write(item[3])
    #    handle_out.write('\n')
    print('Done.')

if __name__ =='__main__':
    main()
