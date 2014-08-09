with open('test') as In:
		Raw=In.read()
Primer=list()
Primer.append(['1','1',Raw[0:20]])
point=20
number=2
fr='F'
for i in range(300):
		if point>=len(Raw):
				Primer.append([Raw[-1:-20],'end','end'])
				break
		point+=350
		Primer.append([Raw[point:(point+20)],number,fr])
		number+=1
		fr='R'
		point+=50
		Primer.append([Raw[point:(point+20)],number,fr])
		number+=1
		fr='F'
for p in Primer:
	print(''.join(['>',str(p[1]),p[2]]))
	print(p[0])





