from glob import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys


def prepare(filelist, data):
    codon = 3
    cut = slice(3, -3)
    head = slice(0, 3)
    start = 'ATG'
    for fasta_file in filelist:
        raw = list(SeqIO.parse(fasta_file, 'fasta'))
        for record in raw:
            length = len(record.seq)
            if length < min_length or length % codon != 0:
                print('{0} in {1} has wrong length.\n'.format(record, 
                                                              fasta_file))
                continue
            name = fasta_file.replace('.fasta', '_' + record)
            if seq(record.seq[head]) != start:
                sequence = record.seq.reverse_complement()
            else:
                sequence = record.seq
            sequence = record.seq[cut]
            data.append(name, sequence)
    return data

def count_gc(filelist,result):
    result = list()
    result.append(['name', 'gc1', 'gc2', 'gc3'])
    ignore = ['ATG', 'TGG']
    gc = 'GC'
    for record in data:
        gc1 = 0
        gc2 = 0
        gc3 = 0
        for i in range(len(data[1]), 3):
            now = data[1][i:i+2]
            if now in ignore:
                continue
            if now[0] in gc:
                gc1 += 1
            if now[1] in gc:
                gc2 += 1
            if now[2] in gc:
                gc3 += 1
        result.append([data[0], gc1, gc2, gc3])
    return result

def main():
    path = sys.argv[1]
    min_length = 300
    filelist = glob(path)
    data = list()
    data.append([name, seq])
    prepare(filelist, data)
    result = count_gc(data)
    with open('result.csv', 'w') as handle:
        for i in result:
            handle.write(','.join(i))
    return 0

if __name__ == '__main__':
    main()
