#!/usr/bin/python3
import sys
import re
from Bio import SeqIO
from Bio.SeqIO.QualityIO import PairedFastaQualIterator as fqout

with open(sys.argv[1],'r') as fna:
    In=fna.read()
with open(sys.argv[2],'r') as qual:
    Out=qual.read()
Rawinfo=re.findall('(?<=\>[0-9]{1,2}-)[0-9a-zA-Z_].*[0-9]{5}',In)
for name in Rawinfo:
    info=name[-15:]
    Out=re.sub(info,name,Out,count=1)
with open(sys.argv[2],'w') as out:
    out.write(Out)

fnaout=list()
qualout=list()
fna=SeqIO.to_dict(SeqIO.parse(sys.argv[1],'fasta'),key_function=lambda rec:rec.id)
qual=SeqIO.to_dict(SeqIO.parse(sys.argv[2],'qual'),key_function=lambda rec:rec.id)
for keys,value in fna.items():
    if keys in qual:
        fnaout.append(value)
        qualout.append(qual[keys])
fnaout.sort(key=lambda x:x.id)
qualout.sort(key=lambda x:x.id)
SeqIO.write(fnaout,sys.argv[1],'fasta')
SeqIO.write(qualout,sys.argv[2],'qual')

fna=sys.argv[1]
qual=sys.argv[2]
fastq=sys.argv[3]
records=fqout(open(fna),open(qual))
SeqIO.write(records,fastq,'fastq')
