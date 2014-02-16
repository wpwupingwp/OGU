from Bio import SeqIO
#File=input("File name:\n")
File="gb"
record=SeqIO.read(File,"genbank")
print(record.id)
exon=[]
for i in record.features:
    if i.type=="gene" and "gene" in i.qualifiers:
        if i.location_operator!="join":
            #Ignore CompoundLocation
            GeneName=str(i.qualifiers["gene"])
            Start=int(i.location.start)
            End=int(i.location.end)
            Strand=int(i.location.strand)
            rec=[GeneName,Start,End,Strand]
            exon.append(rec)
        elif i.location_operator=="join":
            print(i.qualifiers["gene"],i.location)
intron=[]
for i in range(len(exon)-1):
    #Ignore ChongDieGene
    This=exon[i]
    Next=exon[i+1]
    IntronName="-".join([This[0][2:-2],Next[0][2:-2]])
    Start=This[2]+1
    End=Next[1]-1
    print(IntronName,Start,End)


    





