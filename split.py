from Bio import SeqIO
filename=input("filename:\n")
Records=SeqIO.parse(filename,"genbank")
All=[]
for i in Records:
    Taxon=i.features[0].qualifiers["db_xref"][0]
    Accession=i.annotations["Accession"]
    for j in i.features:
        if "gene" in j.qualifiers and "fadR" in j.qualifiers["gene"]:
            Start=j.location.start-400
            End=j.location.end+400
            Sequence=i.seq[Start:End]
            rec=[Taxon,Accession,Sequence]
            All.append(rec)
Out=open(".".join([Accession,"gb"]),"w")
for i in All:
    Out.write(">%s|%s\n%s\n"%(i[0],i[1],i[2]))
print("done.")
