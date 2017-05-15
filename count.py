from Bio import SeqIO
from sys import argv
from glob import glob 
from os import rename

files = glob('*.fa*')
handle = open('count.csv', 'w')
for fastq in files:
    n = 0
    for record in SeqIO.parse(fastq, 'fasta'):
        n += 1
    rename(fastq, '{}-size-{}'.format(fastq, n))
