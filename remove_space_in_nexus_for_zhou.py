#!/usr/bin/python3

from sys import argv
import re


pattern = re.compile(r"(^'?.+'?) +([ATCG-]{20})")
with open(argv[1], 'r') as old, open(argv[1]+'.new', 'w') as new:
    for line in old:
        a = re.search(pattern, line)
        if a is not None:
            old_id = a.group(1)
            new_id = re.sub(r'\W', '_', old_id)
            new_id = '{:<100}'.format(new_id)
            line = line.replace(old_id, new_id)
        new.write(line)
