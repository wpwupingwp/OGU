#!/usr/bin/python3

from Bio import Phylo as p
from sys import argv

try:
    a = p.read(argv[1], 'newick')
except:
    try:
        a = p.read(argv[1], 'nexus')
    except:
        raise Exception('File not recognized.')
b = a.get_terminals()
for i,j in enumerate(b):
    #j.name = 'node{}'.format(i)
    j.name = ''
b_l = len(b)
c = a.get_nonterminals()
c_l = len(c)
p.draw(a)
print('{}\t\t{:.3f}'.format(argv[1], c_l/b_l))

# with open(argv[1], 'r') as raw:
#     data = raw.read()
# node = data.count(',') + 1
# node2 = data.count('(')
# print('{}\t{:.3f}\t{:.3f}'.format(argv[1], c_l/b_l, node2/node))
