#!/usr/bin/python3

import sys

print('Usage: python3 join_fq.py r1.fastq r2.fastq')
with open(sys.argv[1], 'r') as l:
    left = l.read().split(sep='\n')
with open(sys.argv[2], 'r') as r:
    right = r.read().split(sep='\n')
seq = 'NNNNNNNNNN'
qual = 'AAAAAAAAAA'
for i, j in left, right:
    print(i,j)
    


