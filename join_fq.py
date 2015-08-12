#!/usr/bin/python3

import sys

print('Usage: python3 join_fq.py r1.fastq r2.fastq')
with open(sys.argv[1], 'r') as l:
    left = l.read().split(sep='\n')
with open(sys.argv[2], 'r') as r:
    right = r.read().split(sep='\n')
print(left.pop())
print(right.pop())
join_seq = 'NNNNNNNNNN'
join_qual = 'AAAAAAAAAA'
#for i, j,index in enumerate(zip(left, right)):
length = int(len(left)/4)
#every record in fastq file have four lines:
#id\n seq\n id\n qual
handle = open('combine.fastq', 'w')
for index in range(length):
    point = index * 4
    name = left[point]
    seq = ''.join([
        left[point+1],
        join_seq,
        right[point+1]
    ])
    qual_name = left[point+2]
    qual = ''.join([
        left[point+3],
        join_qual,
        right[point+3]
    ])
    handle.write(''.join([
        name,
        '\n',
        seq,
        '\n',
        qual_name,
        '\n',
        qual,
        '\n'
    ]))
