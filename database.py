from Bio.Seq import MutableSeq
import psycopg2
Name=input("Name:\n")
Type=input("Type:\n")
conn=psycopg2.connect(database="data",user="postgres",password="genome",host="127.0.0.1",port="5433")
cur=conn.cursor()
Sql="select Taxon,Organism,Name,Type,Strand,Sequence from main where Name=%s and Type=%s"
cur.execute(Sql,(Name,Type))
Result=cur.fetchall()

All=[]
for i in Result:
    Title="|".join([str(i[0]),i[1],i[2],i[3]])
    Sequence=MutableSeq(i[5])
    if i[4]=="-1":
        Sequence.seq=Sequence.reverse_complement()
    Record=[Title,Sequence]
    All.append(Record)

Fileout=open((".".join(["_".join([Name,Type]),"fasta"])),"w")
for i in All:
    Fileout.write(">%s\n%s\n"%(i[0],i[1]))

cur.close()
conn.close()
