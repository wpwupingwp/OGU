#!/usr/bin/python3

import sys

def main():
    """Expand a given table according to range.
    Input table (CSV format) looks like this:
    >    A,B,C
    It will generate a new table:
    >    D,E 
    where D was expanded from range(B, C) and E is related A.
    """
    print(main.__doc__)
    handle = open('output.csv', 'w')
    with open(sys.argv[1], 'r') as f:
        table = [i.split(sep=',') for i in f.read().split(sep='\n')]
        table.pop(0)
        table.pop(-1)
    pass
    for line in table:
        for i in range(int(line[1]), int(line[2])+1):
            handle.write('{0:06d},{1}\n'.format(i, line[0]))


if __name__ == '__main__':
    main()