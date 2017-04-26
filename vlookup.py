#!/usr/bin/python3

want = set()
with open('./BOP.txt', 'r') as _:
    for line in _:
        want.add(line.split(',')[0].strip())

handle = open('result.csv', 'w')
with open('./DNAbank-v3.csv', 'r') as data:
    for n, line in enumerate(data):
        bop = line.split(',')[0].strip()
        if n == 0:
            handle.write(line)
        if bop in want:
            handle.write(line)
