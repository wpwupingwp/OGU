#!/usr/bin/python3

# Count sequence if barcode in 5' are same in pair-end

from Bio import SeqIO
from glob import glob

F_END = '_1.fastq'
R_END = '_2.fastq'
list_raw = glob('*.fastq')
file_list = {i[:-8] for i in list_raw}
print('Forward,Reverse,Total,Same')
for pair in file_list:
    with open(pair+F_END, 'r') as F, open(pair+R_END, 'r') as R:
        same = 0
        total = 0
        forward = SeqIO.parse(F, 'fastq')
        reverse = SeqIO.parse(R, 'fastq')
        for line in zip(forward, reverse):
            total += 1
            if line[0].seq[:8] == line[1].seq[:8]:
                same += 1
            else:
                pass
                # print(line[0].seq[:8], line[1].seq[:8])
        print(','.join([pair+F_END, pair+R_END, str(total), str(same)]))
