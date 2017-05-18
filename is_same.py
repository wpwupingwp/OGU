#!/usr/bin/python3

from glob import glob
from os import mkdir
from os.path import join as path_join
from shutil import copyfile

mkdir('out')
id_list = {'BOP205013', 'BOP205033', 'BOP205033-NEW', 'BOP205044',
           'BOP205044-NEW'}
file_list = glob('*.fasta')

for fasta in file_list:
    bop = list()
    with open(fasta, 'r') as raw:
        for line in raw:
            if line.startswith('>'):
                bop.append(line[:14])
    is_same = set()
    for item in bop:
        item = item.split('_')
        if item[1] == 'NEW':
            item = '-'.join(item)
        else:
            item = item[0]
        is_same.add(item)
    if len(is_same) == 5:
        copyfile(fasta, path_join('out', fasta))
    else:
        print('\t'.join(is_same))
        print('\t'.join(bop))
