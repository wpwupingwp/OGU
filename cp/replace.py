#!/usr/bin/python3

from glob import glob


info = dict()
with open('info.csv', 'r') as raw:
    for line in raw:
        line = line.split(',')
        line = [i.strip() for i in line]
        info[line[0]] = line[1]

file_list = glob('*.fasta')
for fasta in file_list:
    with open(fasta+'.rename', 'w') as new, open(fasta, 'r') as old:
        for line in old:
            if not line.startswith('>'):
                new.write(line)
            else:
                line = line.split('|')
                newline = list()
                for item in line:
                    try:
                        item = info[item]
                    except:
                        pass
                    newline.append(item)
                new.write('|'.join(newline))
