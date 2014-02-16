from Bio import SeqIO
#FileIn=input("File name:\n")
FileIn="gb"
Record=SeqIO.read(FileIn,"genbank")
Organism=Record.annotations["organism"].replace(" ","_")
Accession=Record.annotations["accessions"][0]
Exon=[]
for i in Record.features:
    if i.type=="gene" and "gene" in i.qualifiers:
        if i.location_operator!="join":
            #Ignore CompoundLocation
            GeneName=str(i.qualifiers["gene"][0])
            Start=int(i.location.start)
            End=int(i.location.end)
            Strand=int(i.location.strand)
            Sequence=str(Record.seq[Start:End])
            rec=[GeneName,Start,End,Strand,Sequence]
            Exon.append(rec)
        elif i.location_operator=="join":
            print(i.qualifiers["gene"],i.location)
Intron=[]
for i in range(len(Exon)-1):
    #Ignore ChongDieGene
    This=Exon[i]
    Next=Exon[i+1]
    IntronName="-".join([This[0],Next[0]])
    Start=This[2]+1
    End=Next[1]-1
    Sequence=str(Record.seq[Start:End])
    rec=[IntronName,Start,End,Sequence]
    Intron.append(rec)
FileOut=open("Out","w")
for i in Exon:
    FileOut.write(">%s|%s|%s\n"%(Organism,i[0],Accession))
    FileOut.write("%s\n"%(i[4]))
for i in Intron:
    FileOut.write(">%s|%s_Intron|%s\n"%(Organism,i[0],Accession))
    FileOut.write("%s\n"%(i[3]))
FileOut.write(">%s|Complete Sequence\n%s"%(Organism,str(Record.seq)))
