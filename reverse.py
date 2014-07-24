#!/usr/bin/python3
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
Records=list(SeqIO.parse(sys.argv[1],"fasta"))
handle=open("reverse_complement.fasta","a")
for item in Records:
    record=SeqRecord(id=item.id,description="",seq=item.seq.reverse_complement())
    SeqIO.write(record,handle,"fasta")
