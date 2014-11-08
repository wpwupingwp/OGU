#!/usr/bin/python3

import sys
from Bio import SeqIO

with open(sys.argv[1],'r') as In:
    Raw=list(SeqIO.parse(In,'gb'))
with open('result.txt','w') as Out:
    for record in Raw:
        Out.write(''.join([record.id,'|',record.description,'\n']))
