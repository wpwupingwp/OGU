import sys
import re
handle1=open(sys.argv[1],"r")
handle2=open(sys.argv[2],"r")
In=str(handle1.read())
Out=str(handle2.read())
handle1.close()
handle2.close()
Rawinfo=re.findall("(?<=\>)[0-9a-zA-Z_].*[0-9]{5}",In)
for name in Rawinfo:
    info=name[-17:]
    Out=re.sub(info,name,Out,count=1)
handle3=open(sys.argv[2],"w")
handle3.write(Out)
