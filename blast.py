from Bio.Blast.Applications import NcbiblastnCommandline as nb
from Bio.Blast import NCBIXML as nx 

cmd=nb(query='unknown.fasta',db='primer',evalue=0.001,outfmt=5,out='O')
stdout,stderr=cmd()
result=nx.parse('O')
