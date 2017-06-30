#!/usr/bin/python3

from Bio import SeqIO
from sys import argv
from timeit import default_timer as timer
import hashlib


def test_format(file_name):
    with open(file_name, 'r') as _:
        line = _.readline()
        if line.startswith('>'):
            return 'fasta'
        elif line.startswith('#'):
            return 'nexus'
        else:
            raise ValueError('Only support fasta and nexus! Quit now.')


def main():
    start = timer()
    file_format = test_format(argv[1])
    raw = SeqIO.parse(argv[1], file_format)
    hash_dict = dict()
    output = list()
    before = 0
    after = 0
    for record in raw:
        before += 1
        h = hashlib.sha512()
        b_text = str(record.seq).encode('utf-8')
        h.update(b_text)
        hash_text = h.hexdigest()
        if hash_text in hash_dict:
            hash_dict[hash_text].append(record)
        else:
            hash_dict[hash_text] = [record, ]
    print('Duplicated sequences:')
    log = open(argv[1]+'.log', 'w')
    for record in hash_dict.values():
        output.append(record[0])
        after += 1
        if len(record) != 1:
            id_list = [i.id for i in record]
            log.write('\t'.join(id_list)+'\n')
    log.write('Before\t{}\tAfter\t{}\n'.format(before, after))
    with open(argv[1]+'.new', 'a') as output_file:
        SeqIO.write(output, output_file, file_format)
    end = timer()
    print('Total {} sequences in {} format.'.format(before, file_format))
    print('{} sequences left in the file {}.new.'.format(after, argv[1]))
    print('Cost {:.3f} seconds.'.format(end-start))


if __name__ == '__main__':
    main()
