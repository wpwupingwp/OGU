#!/usr/bin/python3
import re
In=input("filename:\n")
Adapter=input("adapter:\n")
handle=open(In,"r")
Input=str(handle.read())
out=re.sub(Adapter,"",Input)
handle=open("Out","w")
handle.write(out)
