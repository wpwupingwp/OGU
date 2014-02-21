from Bio import SeqIO
#FileIn=input("File name:\n")
FileIn="gb"
Record=SeqIO.read(FileIn,"genbank")
Organism=Record.annotations["organism"].replace(" ","_")
Accession=Record.annotations["accessions"][0]
Gene=[]
Spacer=[]
rRNA=[]
Exon=[]
Intron=[]
All=[]
for i in Record.features:
    if i.type=="gene" and "gene" in i.qualifiers:
        if i.location_operator!="join":
            Start=int(i.location.start)
            End=int(i.location.end)
            Sequence=str(Record.seq[Start:End])
            Name=str(i.qualifiers["gene"][0])
            Strand=int(i.location.strand)
            rec=[Name,Start,End,Strand,Sequence]
            Gene.append(rec)
        elif i.location_operator=="join":
            Start=int(i.sub_features[0].location.start)
            End=int(i.sub_features[0].location.end)
            Name=str(i.qualifiers["gene"][0])
            Strand=int(i.location.strand)
            Sequence=""
            rec=[Name,Start,End,Strand,Sequence]
            Gene.append(rec)
            Start=int(i.sub_features[1].location.start)
            End=int(i.sub_features[1].location.end)
            Sequence="".join([str(Record.seq[Start:End]),str(Record.seq[Start:End])])
            rec=[Name,Start,End,Strand,Sequence]
            Gene.append(rec)
    elif i.type=="rRNA":
        Start=int(i.location.start)
        End=int(i.location.end)
        Sequence=str(Record.seq[Start:End])
        Name=str(i.qualifiers["product"][0]).replace(" ","_")
        Strand=int(i.location.strand)
        rec=[Name,Start,End,Strand,Sequence]
        rRNA.append(rec)
    elif i.type=="exon":
        Start=int(i.location.start)
        End=int(i.location.end)
        Sequence=str(Record.seq[Start:End])
        Name="_".join([str(i.qualifiers["gene"][0]),"exon",str(i.qualifiers["number"][0])])
        Strand=int(i.location.strand)
        rec=[Name,Start,End,Strand,Sequence]
        Exon.append(rec)
    elif i.type=="intron":
        Start=int(i.location.start)
        End=int(i.location.end)
        Sequence=str(Record.seq[Start:End])
        Strand=int(i.location.strand)
        if "number" in i.qualifiers:
            Name="_".join([str(i.qualifiers["gene"][0]),"intron",str(i.qualifiers["number"][0])])
        else:
            Name="_".join([str(i.qualifiers["gene"][0]),"intron"])
        rec=[Name,Start,End,Strand,Sequence]
        Intron.append(rec)
Gene.sort(key=lambda x:x[1])
for i in range(len(Gene)-1):
    This=Gene[i]
    Next=Gene[i+1]
    Tail=This[2]+1
    Head=Next[1]-1
    Sequence=str(Record.seq[Tail:Head])
    Name="_".join(["-".join([This[0],Next[0]]),"Spacer"])
    Strand=This[3]
    rec=[Name,Tail,Head,Strand,Sequence]
    Spacer.append(rec)
All.extend(Gene)
All.extend(Spacer)
All.extend(Intron)
All.extend(Exon)
All.extend(rRNA)
All.sort(key=lambda x:x[1])
FileOut=open("out","w")
for i in All:
    if i[4]!="":
        FileOut.write(">%s|%s|%s|\n"%(Organism,i[0],Accession,))
        FileOut.write("%s\n"%(i[4]))
FileOut.write(">%s|Complete Sequence\n%s"%(Organism,str(Record.seq)))
