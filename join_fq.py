#!/usr/bin/python3

import sys

print('Usage: python3 join_fq.py r1.fastq r2.fastq')
with open(sys.argv[1], 'r') as l:
    left = l.read().split(sep='\n')
with open(sys.argv[2], 'r') as r:
    right = r.read().split(sep='\n')
print(left.pop())
print(right.pop())
seq = 'NNNNNNNNNN'
qual = 'AAAAAAAAAA'
#for i, j,index in enumerate(zip(left, right)):
length = len(left)/4
print(length)
    
    


