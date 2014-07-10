from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
File=input("file name:\n")
#adapter=input("adapter:\n")
adapter="GACTACGCGTCTAGT"
regex="".join(["^",adapter])
Output=[]
Records=SeqIO.parse(File,"sff")
for record in Records:
    trimed=SeqRecord(id=record.id,description="",seq=record.seq[13:])
    Output.append(trimed)
handle=open(File.replace("sff","fasta"),"w")
SeqIO.write(Output,handle,"fasta")
