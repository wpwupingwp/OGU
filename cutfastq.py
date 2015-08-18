import sys
from Bio import SeqIO

'''
usage: python3 cutfastq.py fastqfile cut_number
'''
maxmium = int(int(sys.argv[2]) / 2)
trim1 = sys.argv[1].replace('.', '-' + sys.argv[2] + '.')
trim2 = sys.argv[1].replace('.', '-2-' + sys.argv[2] + '.')
handle1 = open(trim1, 'a')
handle2 = open(trim2, 'a')
raw = SeqIO.parse(sys.argv[1], 'fastq')
for counts, reads in enumerate(raw):
    if counts < maxmium:
        SeqIO.write(reads, handle1, 'fastq')
    else:
        SeqIO.write(reads, handle2, 'fastq')
