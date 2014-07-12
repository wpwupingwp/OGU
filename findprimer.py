#draft
from Bio import SeqIO
import re
Sequence=input("sequence name:\n")
Primer=input("primer list name:\n")
Records=list(SeqIO.parse(Sequence,"fasta"))
Primers=list(SeqIO.parse(Primer,"fasta"))
Output=[]
for record in Records:
    for primer in Primers:
        if re.search(str(primer.seq),str(record.seq)):
            Output.append(record,primer)
            break
for item in Output:
    SeqIO.write(item[0],item[1],"fasta")




