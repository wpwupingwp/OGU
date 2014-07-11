from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
File=input("file name:\n")
Output=[]
Records=SeqIO.parse(File,"fasta")
for record in Records:
    Reverse=SeqRecord(id=record.id,description="",seq=record.seq)
    Output.append(Reverse)
handle=open("reverse.fasta","w")
SeqIO.write(Output,handle,"fasta")
