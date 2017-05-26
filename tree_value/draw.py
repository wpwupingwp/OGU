#!/usr/bin/python3

a = p.read('./mega_bp.nwk', 'newick')
b = a.get_terminals()
b_l = len(b)
c = a.get_nonterminals()
c_l = len(c)
print('{:.3f}'.format(c_l/b_l))
