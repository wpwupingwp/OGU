#!/usr/bin/python3

a = dict()
with open('./Takhtajan.txt', 'r') as _:
    for line in _:
        line = line.strip().split('\t')
        a[line[0]] = line

with open('./APG4.txt', 'r') as b:
    for line in b:
        head = line.split('\t')[0]
        if head in a:
            a[head].append(line.strip())

with open('./Thorne.txt', 'r') as c:
    for line in c:
        head = line.split('\t')[0]
        if head in a:
            a[head].append(line.strip())

with open('out.txt', 'w') as output:
    for i in a.values():
        print(i)
        output.write('{}\n'.format('\t'.join(i)))
