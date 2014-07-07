from Bio import SeqIO
File=input("File name:\n")
a=SeqIO.parse(File,"fasta")
List=[]
for i in a:
    List.append([i.id,len(i.seq)])
d=0
e=0
for i in List:
    print(i[0],i[1])
    d+=1
    e+=i[1]
print(len(List),e/float(len(List)))

