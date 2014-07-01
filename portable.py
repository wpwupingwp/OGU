from Bio import SeqIO,BiopythonDeprecationWarning
from Bio.Seq import MutableSeq
import sqlite3,warnings,datetime
warnings.simplefilter("ignore",BiopythonDeprecationWarning)

def parser():
    Taxon=int(Record.features[0].qualifiers["db_xref"][0][6:])
    Organism=Record.annotations["organism"]
    Accession=Record.annotations["accessions"][0]
    Gene=[]
    All=[]
    for i in Record.features:
        if i.type=="gene" and "gene" in i.qualifiers:
            if i.location_operator!="join":
                Type="gene"
                Start=int(i.location.start)
                End=int(i.location.end)
                Sequence=str(Record.seq[Start:End])
                Name=str(i.qualifiers["gene"][0])
                Strand=str(i.location.strand)
                rec=[Taxon,Organism,Accession,Name,Type,Start,End,Strand,Sequence,Date]
                Gene.append(rec)
            elif i.location_operator=="join":
                Type="gene"
                Start=int(i.sub_features[0].location.start)
                End=int(i.sub_features[0].location.end)
                Name=str(i.qualifiers["gene"][0])
                Strand=str(i.location.strand)
                Sequence=""
                rec=[Taxon,Organism,Accession,Name,Type,Start,End,Strand,Sequence,Date]
                Gene.append(rec)
                Start=int(i.sub_features[1].location.start)
                End=int(i.sub_features[1].location.end)
                Sequence="".join([str(Record.seq[Start:End]),str(Record.seq[Start:End])])
                rec=[Taxon,Organism,Accession,Name,Type,Start,End,Strand,Sequence,Date]
                Gene.append(rec)
        elif i.type=="rRNA":
            Type="rRNA"
            Start=int(i.location.start)
            End=int(i.location.end)
            Sequence=str(Record.seq[Start:End])
            Name=str(i.qualifiers["product"][0]).replace(" ","_")
            Strand=str(i.location.strand)
            rec=[Taxon,Organism,Accession,Name,Type,Start,End,Strand,Sequence,Date]
            All.append(rec)
        elif i.type=="exon" and "gene" in i.qualifiers :
            Type="exon"
            Start=int(i.location.start)
            End=int(i.location.end)
            Sequence=str(Record.seq[Start:End])
            if "number" in i.qualifiers:
                Name="_".join([str(i.qualifiers["gene"][0]),"exon",str(i.qualifiers["number"][0])])
            else:
                Name="_".join([str(i.qualifiers["gene"][0]),"exon"])
            Strand=int(i.location.strand)
            rec=[Taxon,Organism,Accession,Name,Type,Start,End,Strand,Sequence,Date]
            All.append(rec)
        elif i.type=="intron" and "gene" in i.qualifiers:
            Type="intron"
            Start=int(i.location.start)
            End=int(i.location.end)
            Sequence=str(Record.seq[Start:End])
            Strand=str(i.location.strand)
            if "number" in i.qualifiers:
                Name="_".join([str(i.qualifiers["gene"][0]),"intron",str(i.qualifiers["number"][0])])
            else:
                Name="_".join([str(i.qualifiers["gene"][0]),"intron"])
            rec=[Taxon,Organism,Accession,Name,Type,Start,End,Strand,Sequence,Date]
            All.append(rec)
    Gene.sort(key=lambda x:x[5])
    for i in range(len(Gene)-1):
        Type="spacer"
        This=Gene[i]
        Next=Gene[i+1]
        Tail=This[6]+1
        Head=Next[5]-1
        Sequence=str(Record.seq[Tail:Head])
        Name="_".join(["-".join([This[3],Next[3]]),"Spacer"])
        Strand=0
        rec=[Taxon,Organism,Accession,Name,Type,Start,End,Strand,Sequence,Date]
        All.append(rec)
    All.extend(Gene)
    Database.extend(All)
    return 


def database():
    con=sqlite3.connect("./db")
    cur=con.cursor()
    cur.execute("create table if not exists main (Taxon int,Organism text,Accession text,Name text,Type text,Head int,Tail int, Strand text,Sequence text,Date text,ID integer PRIMARY KEY);")
    for row in Database:
        if row[9]!="":
            cur.execute("insert into main (Taxon,Organism,Accession,Name,Type,Head,Tail,Strand,Sequence,Date) values (?,?,?,?,?,?,?,?,?,?);",(row[0],row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8],row[9]))
    con.commit()
    cur.close()
    con.close()
    print("Done.\n")
    return
    
def Query():
    Querytype=input("1.Specific fragment\n2.Specific Organism\n3.Specific gene\n4.All\n")
    if Querytype in ["1","2","3","4"]:
        RunQuery(Querytype)
    else:
        print("Input error!\n")
    return

def RunQuery(Querytype):
    con=sqlite3.connect("./db")
    cur=con.cursor()
    if Querytype=="1":
        Organism=input("Organism:\n")
        Gene=input("Gene:\n")
        Type=input("Fragment type(gene,rRNA,exon,intron,spacer):\n")
        cur.execute("select Taxon,Organism,Name,Type,Strand,Sequence from main where Name like ? and Type=? and Organism like ?",('%'+Gene+'%',Type,Organism))
        Result=cur.fetchall()
    elif Querytype=="2":
        Organism=input("Organism:\n")
        cur.execute("select Taxon,Organism,Name,Type,Strand,Sequence,Head from main where Organism=? order by Head",(Organism,))
        Result=cur.fetchall()
    elif Querytype=="3":
        Gene=input("Gene:\n")
        Type=input("Fragment type(gene,rRNA,exon,intron,spacer):\n")
        cur.execute("select Taxon,Organism,Name,Type,Strand,Sequence from main where Name like ? and Type=? order by Taxon",('%'+Gene+'%',Type))
        Result=cur.fetchall()
    elif Querytype=="4":
        cur.execute("select Taxon,Organism,Name,Type,Strand,Sequence,Head from main order by Taxon")
        Result=cur.fetchall()

    All=[]
    for i in Result:
        Title="|".join([str(i[0]),i[1],i[2],i[3]])
        Sequence=MutableSeq(i[5])
        if i[4]=="-1":
            Sequence.seq=Sequence.reverse_complement()
        Record=[Title,Sequence]
        All.append(Record)

    Output=input("Enter output filename:\n")
    Fileout=open(".".join([Output,"fasta"]),"w")
    for i in All:
        Fileout.write(">%s\n%s\n"%(i[0],i[1]))
    cur.close()
    con.close()
    print("Done.\n")
    return 

#Main program 
Option=input("Select:\n1.Add data\n2.Query\n")
Database=[]
Date=str(datetime.datetime.now())[:19].replace(" ","-")
if Option=="1":
   FileIn=input("Genbank format filename:\n")
   Records=SeqIO.parse(FileIn,"genbank")
   for Record in Records:
       parser()
   database()
elif Option=="2":
    Query()
else:
    print("Input error!\n")

