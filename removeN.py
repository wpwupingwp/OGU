from Bio import SeqIO
print("This little program will remove the sequence which have too many N(one percent).\n")
In=input("Filename:\n")
handle2=open("output","w")
Input=list(SeqIO.parse(In,"fasta"))
Output=[]
for record in Input:
    N=0
    Len=len(record.seq)
    for letter in str(record.seq):
        if letter in ["M","R","S","W","Y","K","N"]:
            N+=1
    if N/float(Len)>=0.01:
        print(record,N,Len)
    else:
        Output.append(record)
SeqIO.write(Output,"output","fasta")
print(len(Input)-len(Output),"sequences removed.\n")

