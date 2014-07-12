import re
In=input("filename")
handle=open(In,"r")
Input=str(handle.read())
out=re.sub("GACTACGCGTCTAGT","",Input)
handle=open("Out","w")
handle.write(out)
