from Bio import SeqIO
#FileIn=input("File name:\n")
FileIn="gb"
Record=SeqIO.read(FileIn,"genbank")
Organism=Record.annotations["organism"].replace(" ","_")
Accession=Record.annotations["accessions"][0]
Exon=[]
for i in Record.features:
    if i.type=="gene" and "gene" in i.qualifiers:
        #Ignore genes without name
        if i.location_operator!="join":
            Join=0
            Start=[int(i.location.start)]
            End=[int(i.location.end)]
            Start.append("-1")
            End.append("-1")
            Sequence=str(Record.seq[Start[0]:End[0]])
            GeneName=str(i.qualifiers["gene"][0])
            Strand=int(i.location.strand)
            rec=[GeneName,Join,Start,End,Strand,Sequence]
            Exon.append(rec)
        elif i.location_operator=="join":
            Join=1
            Start=[int(i.sub_features[0].location.start)]
            End=[int(i.sub_features[0].location.end)]
            Start.append(int(i.sub_features[1].location.start))
            End.append(int(i.sub_features[1].location.end))
            GeneName=str(i.qualifiers["gene"][0])
            Strand=int(i.location.strand)
            Sequence="".join([str(Record.seq[Start[0]:End[0]]),str(Record.seq[Start[1]:End[1]])])
            rec=[GeneName,Join,Start,End,Strand,Sequence]
            Exon.append(rec)
Intron=[]
for i in range(len(Exon)-1):
    #Ignore Overlap gene 
    This=Exon[i]
    Next=Exon[i+1]
    IntronName="-".join([This[0],Next[0]])
    Start=This[2][0]+1
    End=Next[3][0]-1
    Sequence=str(Record.seq[Start:End])
    rec=[IntronName,Start,End,Sequence]
    Intron.append(rec)
FileOut=open("out","w")
for i in Exon:
#    print(i)
    FileOut.write(">%s|%s|%s|\n"%(Organism,i[0],Accession,))
    FileOut.write("%s\n"%(i[5]))
for i in Intron:
    FileOut.write(">%s|%s_Intron|%s\n"%(Organism,i[0],Accession))
    FileOut.write("%s\n"%(i[3]))
FileOut.write(">%s|Complete Sequence\n%s"%(Organism,str(Record.seq)))
