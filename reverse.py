from Bio import SeqRecord,SeqIO
import sys
Records=list(SeqIO.parse(sys.argv[1],"fasta"))
handle=open("trim.fasta","a")
for item in Records:
    record=SeqRecord(id=item.id,description="",seq=item.seq.reverse_complement())
    SeqIO.write(record,"fasta")
