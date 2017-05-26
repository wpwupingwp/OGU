#!/usr/bin/python3

from Bio import Phylo as p
from sys import argv

a = p.read(argv[1], 'newick')
b = a.get_terminals()
b_l = len(b)
c = a.get_nonterminals()
c_l = len(c)
p.draw(a)
print('{:.3f}'.format(c_l/b_l))
