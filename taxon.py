#!/usr/bin/python3
import sqlite3

def Create():
    Id=dict()
    Rank=dict()
    Name=dict()
    Species=[]
    NotSpecies=[]
    Data=[]
    with open('./test/nodes','r') as In:
        Raw=list(In.readlines())
        for record in Raw:
            add=record.replace('\n','').split(sep=' ')
            Id[add[0]]=add[1]
            Rank[add[0]]=add[2]
            if add[2]=='species':
                Species.append(add[0])
            else:
                NotSpecies.append(add[0])
    with open('./test/name','r') as In:
        Raw=list(In.readlines())
        for record in Raw:
            add=record.replace('\n','').split(sep='|')
            if add[0] not in Name or add[3]=='scientific name':
                Name[add[0]]=add[1]
    for species in Species:
        record=[species,]
        while Id[species]!='1' :
    #        if Rank[Id[species]] in ['species','genus','family','order','class','phylum','kingdom'] :
            record.append(Id[species])
            species=Id[species]
        if '33090' in record:
            record.pop()
            record.pop()
            Data.append(record[::-1])
    return


def database():
    con=sqlite3.connect("./test/DB")
    cur=con.cursor()
    cur.execute("create table if not exists name_taxonid (ID integer PRIMARY KEY,Name text);")
    for row in Name:
        cur.execute("insert into name_taxonid (ID,Name) values (?,?);",(row[0],row[1]))
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
        Type=input("Fragment type(gene,cds,rRNA,tRNA,exon,intron,spacer):\n")
        cur.execute("select Taxon,Organism,Name,Type,Strand,Sequence from main where Name like ? and Type=? and Organism=?",('%'+Gene+'%',Type,Organism))
        Result=cur.fetchall()
    elif Querytype=="2":
        Organism=input("Organism:\n")
        Type=input("Fragment type(gene,cds,rRNA,tRNA,exon,intron,spacer,whole,fragments):\n")
        if Type=="fragments":
            cur.execute("select Taxon,Organism,Name,Type,Strand,Sequence,Head from main where Organism=?  order by Head",(Organism,))
        else:
            cur.execute("select Taxon,Organism,Name,Type,Strand,Sequence,Head from main where Organism=? and Type=? order by Head",(Organism,Type))
        Result=cur.fetchall()
    elif Querytype=="3":
        Gene=input("Gene:\n")
        Type=input("Fragment type(gene,cds,rRNA,tRNA,exon,intron,spacer):\n")
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
    Fileout.close()
    print("Done.\n")
    return 

