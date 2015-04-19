from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys

in_name=sys.argv[1]
out_name=sys.argv[1].replace('.','_r.')
raw=list(SeqIO.parse(in_name,'fasta'))
out=list()
for n in range(len(raw)):
    sequence=raw[n].seq
    new_id=''.join(['scaffold',str(n+1),'.1|size',str(len(sequence))])
    to_add=SeqRecord(id=new_id,description='',seq=sequence)
    out.append(to_add)
SeqIO.write(out,out_name,'fasta')

