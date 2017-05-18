#!/usr/bin/python3

import re
from glob import glob
from os import rename

file_list = glob('*.*')
pattern = re.compile(r'size[-_](\d+)')
for fasta in file_list:
    total = 0
    with open(fasta, 'r') as raw:
        for line in raw:
            if line.startswith('>'):
                match = pattern.search(line)
                if match is not None:
                    total += int(match.group(1))
    rename(fasta, '{}-sum-{}.fasta'.format(fasta, total))
