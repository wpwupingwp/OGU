#!/usr/bin/python3

from sys import argv
import re


pattern = re.compile(r'([A-Za-z]) ([A-Za-z])')
with open(argv[1], 'r') as old, open(argv[1]+'.new', 'w') as new:
    for line in old:
        if line.startswith("'"):
            new_line = re.sub(pattern, r'\1_\2', line)
            new.write(new_line)
        else:
            new.write(line)
