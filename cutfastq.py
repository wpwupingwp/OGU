import sys
from Bio import SeqIO

maxmium=600000
#after that, it's very slow (several days)

raw=list(SeqIO.parse(sys.argv[1],'fastq'))
trim=raw[:maxmium]
SeqIO.write(trim,sys.argv[1].replace('.','_cut.'),'fastq')

