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
    duplicate = dict()
    output = list()
    for record in raw:
        h = hashlib.sha512()
        b_text = str(record.seq).encode('utf-8')
        h.update(b_text)
        hash_text = h.hexdigest()
        if hash_text in hash_dict:
            if hash_text in duplicate:
                duplicate[hash_text].append(record.id)
            else:
                duplicate[hash_text] = [record.id, ]
        else:
            output.append(record)
    with open(argv[1]+'.new', 'w') as output_file:
        SeqIO.write(output, output_file, file_format)


if __name__ == '__main__':
    main()
