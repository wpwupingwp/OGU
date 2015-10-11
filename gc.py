from glob import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys


def prepare(filelist, data):
    codon = 3
    min_length = 300
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
            name = record.id
            if str(record.seq[head]) != start:
                sequence = record.seq.reverse_complement()
            else:
                sequence = record.seq
                print(sequence)
            sequence = str(record.seq[cut].upper())
            data.append([fasta_file, name, sequence])
    return data

def count_gc(data):
    result = list()
    result.append(['file', 'name', 'gc1', 'gc2', 'gc3'])
    ignore = ['ATG', 'TGG']
    gc = 'GC'
    for record in data:
        gc1 = 0
        gc2 = 0
        gc3 = 0
        for i in range(0, len(record[2]), 3):
            now = record[2][i:i+3]
            if now in ignore:
                continue
            if now[0] in gc:
                gc1 += 1
            if now[1] in gc:
                gc2 += 1
            if now[2] in gc:
                gc3 += 1
        result.append([record[0], record[1], gc1, gc2, gc3])
    return result

def main():
    path = sys.argv[1]
    filelist = glob(path)
    data = list()
    prepare(filelist, data)
    result = count_gc(data)
   # print(data)
    with open('result.csv', 'w') as handle:
        for i in result:
            handle.write('{0},{1},{2},{3},{4}\n'.format(*i))
    return 0

if __name__ == '__main__':
    main()
