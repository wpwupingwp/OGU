#!/usr/bin/python3
# -*- coding:utf-8 -*-
a=input("Name:\n")
b=input("Enter sequence:\n")
b=b.upper()
C=b.count("C")
T=b.count("T")
Tm=4*2*C+2*2*T
Start=b.find("AUG")
print("\n")
print("The species is",a)
print("Its sequence is:\n",b.upper(),"\n",sep="")
print("Its length is",b.__len__())
print("Its Tm value may be",Tm)
print("Its start codon is in",Start)
print("Its mRNA is",b.replace("T","U"))
print(b.split("A"),"\n",b.split("T"),"\n",b.split("C"),"\n",b.split("G"),"\n")
print("\n".join([a,b]))




