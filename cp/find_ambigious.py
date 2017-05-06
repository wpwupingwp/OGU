#!/usr/bin/python3


from Bio import SeqIO
from glob import glob
# try type annotation
from typing import *


# ignore 'N'
letters = set('RYKMSWBDHV')
file_list = glob('*.fasta')
handle = open('ambigious_base.txt', 'w')
handle2 = open('duplicate_sequence.txt', 'w')

for fasta in file_list:
    handle.write(fasta+'\n')
    handle2.write(fasta+'\n')
    id_list:Set[str] = set()
    for index, sequences in enumerate(SeqIO.parse(fasta, 'fasta')):
        if sequences.id in id_list:
            handle2.write('{}. {}\n'.format(index, sequences.id))
        id_list.add(sequences.id)
        for n, letter in enumerate(sequences.seq):
            if letter in letters:
                handle.write(sequences.id+'\n')
                handle.write('Position:\t{0}\tBase:\t{1}\n'.format(
                    n+1, letter))
                break
    handle.write('\n')
    handle2.write('\n')
handle.close()
handle2.close()
