#!/usr/bin/python3
import sys
from Bio import SeqIO
print("python3",sys.argv[0],"sff_file","fasta_file")
SeqIO.convert(sys.argv[1],"sff",sys.argv[2],"fastq")
