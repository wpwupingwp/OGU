import re

Sum=list()
with open('listall.txt','r') as In:
    List=In.read().split('\n')
    List.pop()
with open('log','r') as In:
    Raw=In.read()
for item in List:
    number=len(re.findall(item,Raw))
    Sum.append([item,number])
for item in Sum:
    print(item[0],item[1])
