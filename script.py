from Bio import SeqIO
#File=input("File name:\n")
File="gb"
record=SeqIO.read(File,"genbank")
print(record.id)
for i in record.features:
    if i.type=="gene":
        if "gene" in i.qualifiers:
            print(i.qualifiers["gene"],i.location)
    elif i.type=="intron":
        print(i.location,"intron")
