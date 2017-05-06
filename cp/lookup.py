#!/usr/bin/python3

import argparse


def reformat(want):
    if want[0].startswith('bop'):
        want = [i.upper() for i in want]
    elif want[0].startswith('BOP'):
        pass
    else:
        want = ['BOP'+i for i in want]
    return {i: 0 for i in want}


def query():
    want = list()
    if arg.query_file is not None:
        with open(arg.query_file, 'r', encoding='utf-8') as raw:
            for item in raw.read().split(','):
                want.append(item.strip())
    if arg.query_list:
        raw = input('Input bop list, seperate by English comma:\n')
        for item in raw.split(','):
            want.append(item.strip())
    want = reformat(want)

    handle = open(arg.output, 'w', encoding='utf-8')
    with open(DB, 'r', encoding='utf-8') as data:
        for n, line in enumerate(data):
            bop = line.split(',')[0].strip()
            if n == 0:
                print(line)
                handle.write(line)
            if bop in want:
                print(line)
                handle.write(line)
                want[bop] += 1
    handle.close()
    for name, n in want.items():
        if n == 0:
            print('{0} not found!'.format(name))
        elif n > 1:
            print('Found duplicate records for {0}!'.format(name))
        else:
            pass


def main():
    global arg
    arg = argparse.ArgumentParser()
    arg.add_argument('-q', dest='query_list', action='store_true',
                     help='query id, seperate by English comma')
    arg.add_argument('-l', dest='query_file',
                     help='query id list file, csv format')
    arg.add_argument('-o', '--output', default='result.csv',
                     help='result file')
    arg.print_help()
    arg = arg.parse_args()
    print(arg)
    global DB
    DB = 'DNABank-v3.csv'
    # start here
    if arg.query_file is None and arg.query_list is None:
        raise Exception('Wrong input!')
    query()
    print('Done. Output was written into {0}'.format(arg.output))
    # end


if __name__ == '__main__':
    main()
