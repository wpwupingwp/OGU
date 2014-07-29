#!/usr/bin/python3
import re
import sys

A1='(gactactcgcgtcgt'
A2='gactacgtacacact'
A3='gactactatacgagt'
A4='gactactacgtctct'
A5='gactacgcgtctagt'
A6='gactacgtagatcgt'
A7=A3
A81='gactacgtactgtgt'
A82='gactacacgacgact'
A91=A5
A92='gactacacgtagtat'
A10=A82
A11='gactacgagtagact'
A12='gactacgacacgtat)'
exp='|'.join([A1,A2,A3,A4,A5,A6,A81,A82,A92,A11,A12])
with open(sys.argv[1],"r") as In:
    Input=str(In.read())
out=re.sub(exp,"",Input)
File=(sys.argv[1]).replace('.','-trim.')
handle=open(File,"w")
handle.write(out)
