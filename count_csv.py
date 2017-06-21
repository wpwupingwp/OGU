#!/usr/bin/python

from sys import argv

raw = list()
with open(argv[1], 'r') as input_file:
    for line in input_file:
        line = line.strip()
        raw.append(line.split(','))
head = {i[0]: list() for i in raw}
print(head)
for column in range(len(raw[0])):
    if column == 0:
        continue
    for l

