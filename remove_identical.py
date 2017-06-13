#!/usr/bin/python3

from Bio import SeqIO
from sys import argv
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
    file_format = test_format(argv[1])
    raw = SeqIO.parse(argv[1], file_format)
    hash_dict = dict()
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
    for record in hash_dict.values():
        if len(record) == 1:
            after += 1
            with open(argv[1]+'.new', 'a') as output:
                SeqIO.write(record, output, file_format)
        else:
            id_list = [i.id for i in record]
            print('Duplicated sequences:')
            print('\t'.join(id_list))
    print('Total {} sequences in {} format.'.format(before, file_format))
    print('{} sequences left in the file {}.new.'.format(after, argv[1]))


if __name__ == '__main__':
    main()
