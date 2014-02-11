#!/usr/bin/python3
# This little program is used to accept inputs and print a FASTA-like record.
a=input("Name:\n")
b=input("Sequence:\n")
print(">",a,"\n",b,sep="",end="\n\n")
print(">"+a+"\n"+b+"\n")


