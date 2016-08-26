#!/usr/bin/python3

import sys
import time

def main():
    '''This file will do something VLOOKUP cannot. It will lookup values in a
    given range defined by two column.
    Table 1 looks like this:
    A, B
    Table 2 looks like this:
    C, D, E
    It will fill B with C if A smaller than E and bigger than D.

    Usage:
    python3 vlookup_assistant.py table_1_file table_2_file
    '''
    print(main.__doc__)
    file_1 = sys.argv[1]
    file_2 = sys.argv[2]
    table_2 = list()
    with open(file_1, 'r') as f:
        table_1 = [[i, ''] for i in f.read().split(sep='\n')]
        table_1.pop()
    with open(file_2, 'r') as f:
        table_2 = [i.split(sep=',') for i in f.read().split(sep='\n')]
        table_2.pop()
    pass

if __name__ == '__main__':
    main()