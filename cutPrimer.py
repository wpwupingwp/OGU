with open('in.fasta') as In:
    Raw=In.read()
Primer=list()
Primer.append([Raw[0:20],1,'F'])
Primer.append([Raw[450:470],2,'F'])
point=450
number=1
fr='R'
for i in range(300):
    if point>=len(Raw):
        break
    point+=30
    p=Raw[point:(point+20)]
    add=[p,number,fr]
    Primer.append(add)
    number+=2
    fr='F'
    point+=420
    p=Raw[point:(point+20)]
    add=[p,number,fr]
    Primer.append(add)
    number-=1
    fr='R'

handle=open('primer.fasta','a')
for p in Primer:
    handle.write(''.join(['>',str(p[1]),p[2]]))
    handle.write('\n')
    handle.write(p[0])
    handle.write('\n')
