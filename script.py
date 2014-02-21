from Bio import SeqIO
#FileIn=input("File name:\n")
FileIn="gb"
Record=SeqIO.read(FileIn,"genbank")
Organism=Record.annotations["organism"].replace(" ","_")
Accession=Record.annotations["accessions"][0]
Gene=[]
Spacer=[]
rRNA=[]
for i in Record.features:
    if i.type=="gene" and "gene" in i.qualifiers:
        #Ignore genes without name
        if i.location_operator!="join":
            Start=int(i.location.start)
            End=int(i.location.end)
            Sequence=str(Record.seq[Start:End])
            GeneName=str(i.qualifiers["gene"][0])
            Strand=int(i.location.strand)
            rec=[GeneName,Start,End,Strand,Sequence]
            Gene.append(rec)
        elif i.location_operator=="join":
            Start=int(i.sub_features[0].location.start)
            End=int(i.sub_features[0].location.end)
            GeneName=str(i.qualifiers["gene"][0])
            Strand=int(i.location.strand)
            Sequence=""
            rec=[GeneName,Start,End,Strand,Sequence]
            Gene.append(rec)
            Start=int(i.sub_features[1].location.start)
            End=int(i.sub_features[1].location.end)
            Sequence="".join([str(Record.seq[Start:End]),str(Record.seq[Start:End])])
            rec=[GeneName,Start,End,Strand,Sequence]
            Gene.append(rec)
    elif i.type=="rRNA":
        Start=int(i.location.start)
        End=int(i.location.end)
        Sequence=str(Record.seq[Start:End])
        rRNAName=str(i.qualifiers["product"][0]).replace(" ","_")
#        Strand=int(i.location.strand)
        rec=[rRNAName,Start,End,Sequence]
        rRNA.append(rec)
        print(rec)
Gene.sort(key=lambda x:x[1])
for i in range(len(Gene)-1):
    This=Gene[i]
    Next=Gene[i+1]
    Tail=This[2]+1
    Head=Next[1]-1
    Sequence=str(Record.seq[Tail:Head])
    SpacerName="-".join([This[0],Next[0]])
    rec=[SpacerName,Tail,Head,Sequence]
    Spacer.append(rec)
FileOut=open("out","w")
for i in Gene:
    if i[4]!="":
        FileOut.write(">%s|%s|%s|\n"%(Organism,i[0],Accession,))
        FileOut.write("%s\n"%(i[4]))
for i in Spacer:
    if i[3]!="":
        FileOut.write(">%s|%s_Spacer|%s\n"%(Organism,i[0],Accession))
        FileOut.write("%s\n"%(i[3]))
for i in rRNA:
    FileOut.write(">%s|%s|%s"%(Organism,i[0],Accession))
    FileOut.write("%s\n"%(i[3]))
FileOut.write(">%s|Complete Sequence\n%s"%(Organism,str(Record.seq)))
