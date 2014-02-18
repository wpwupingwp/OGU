from Bio import SeqIO
#FileIn=input("File name:\n")
FileIn="gb"
Record=SeqIO.read(FileIn,"genbank")
Organism=Record.annotations["organism"].replace(" ","_")
Accession=Record.annotations["accessions"][0]
Exon=[]
Intron=[]
for i in Record.features:
    if i.type=="gene" and "gene" in i.qualifiers:
        #Ignore genes without name
        if i.location_operator!="join":
            Join=0
            Start=int(i.location.start)
            End=int(i.location.end)
            Sequence=str(Record.seq[Start:End])
            GeneName=str(i.qualifiers["gene"][0])
            Strand=int(i.location.strand)
            rec=[GeneName,Join,Start,End,Strand,Sequence]
            Exon.append(rec)
        elif i.location_operator=="join":
            Join=1
            Start=int(i.sub_features[0].location.start)
            End=int(i.sub_features[0].location.end)
            GeneName=str(i.qualifiers["gene"][0])
            Strand=int(i.location.strand)
            Sequence=""
            rec=[GeneName,Join,Start,End,Strand,Sequence]
            Exon.append(rec)
            Start=int(i.sub_features[1].location.start)
            End=int(i.sub_features[1].location.end)
            Sequence="".join([str(Record.seq[Start:End]),str(Record.seq[Start:End])])
            rec=[GeneName,Join,Start,End,Strand,Sequence]
            Exon.append(rec)
Exon.sort(key=lambda x:x[2])
for i in range(len(Exon)-1):
    #Ignore Overlap gene 
    This=Exon[i]
    Next=Exon[i+1]
    Tail=This[3]+1
    Head=Next[2]-1
    Sequence=str(Record.seq[Tail:Head])
    IntronName="-".join([This[0],Next[0]])
    rec=[IntronName,Tail,Head,Sequence]
    Intron.append(rec)
FileOut=open("out","w")
for i in Exon:
#    print(i)
    FileOut.write(">%s|%s|%s|\n"%(Organism,i[0],Accession,))
    if i[5]!="":
        FileOut.write("%s\n"%(i[5]))
for i in Intron:
    FileOut.write(">%s|%s_Intron|%s\n"%(Organism,i[0],Accession))
    if i[3]!="":
        FileOut.write("%s\n"%(i[3]))
FileOut.write(">%s|Complete Sequence\n%s"%(Organism,str(Record.seq)))
