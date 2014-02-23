from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import psycopg2
Name=input("Name:\n")
Type=input("Type:\n")
conn=psycopg2.connect(database="data",user="postgres",password="genome",host="127.0.0.1",port="5433")
cur=conn.cursor()
Sql="select Taxon,Organism,Name,Type,Strand,Sequence from main where Name=%s and Type=%s"
cur.execute(Sql,(Name,Type))
Result=cur.fetchall()

Fileout=(".".join(["_".join([Name,Type]),"fasta"]))
All=[]
for i in Result:
    Title="|".join([i[1],i[2],i[3]])
    Id=str(i[0])
    Sequence=Seq(i[5])
    Record=SeqRecord(Sequence,id=Id,description=Title)
    if i[4]=="-1":
        Record.Seq=Record.reverse_complement
    All.append(Record)
SeqIO.write(All,Fileout,"fasta")

cur.close()
conn.close()
