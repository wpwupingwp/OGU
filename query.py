from Bio.Seq import MutableSeq
import sqlite3

def query():
    con=sqlite3.connect("./db")
    cur=con.cursor()
    Name=input("Name:\n")
    Type=input("Type:\n")
    cur.execute("select Taxon,Organism,Name,Type,Strand,Sequence from main where Name='%s' and Type='%s'" %(Name,Type))
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
    con.close()
    return 

query()
