#!/usr/bin/python

from sys import argv

sample_info = dict()
with open(argv[1], 'r') as input_file:
    is_head = True
    for line in input_file:
        line = line.strip()
        line = line.split(',')
        if is_head:
            head = line
            is_head = False
        else:
            for column in range(len(head)):
                item = line[column]
                if item not in sample_info:
                    sample_info[item] = [[i, 0] for i in head]
                sample_info[item][column][1] = 1
sample_info.pop('')
print('Sample', '\t'.join(head))
for i in sample_info.keys():
    line = [j[1] for j in sample_info[i]]
    print(i, line)
