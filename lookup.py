#!/usr/bin/python3

import argparse
from functools import wraps
from timeit import default_timer as timer


def print_time(function):
    @wraps(function)
    def wrapper(*args, **kargs):
        start = timer()
        result = function(*args, **kargs)
        end = timer()
        print('The function {0} Cost {1:3f}s.\n'.format(
            function.__name__, end-start))
        return result
    return wrapper


@print_time
def reformat(want):
    if want[0].startswith('bop'):
        want = [i.upper() for i in want]
    elif want[0].startswith('BOP'):
        pass
    else:
        want = ['BOP'+i for i in want]
    return set(want)


@print_time
def query():
    want = list()
    if arg.query_file is not None:
        with open(arg.query_list, 'r') as raw:
            for item in raw.read().split(','):
                want.append(item.strip())
    if arg.query_list is not None:
        for item in arg.query_list.split(','):
            want.append(item.strip())
    want = reformat(want)
    print(want)

    handle = open(arg.output, 'w')
    with open(DB, 'r') as data:
        for n, line in enumerate(data):
            bop = line.split(',')[0].strip()
            if n == 0:
                handle.write(line)
            if bop in want:
                handle.write(line)
    handle.close()


def main():
    global arg
    arg = argparse.ArgumentParser()
    arg.add_argument('-q', dest='query_list',
                     help='query id, seperate by English comma')
    arg.add_argument('-l', dest='query_file',
                     help='query id list file, csv format')
    arg.add_argument('-o', '--output', default='result.csv',
                     help='result file')
    arg.print_help()
    arg = arg.parse_args()
    global DB
    DB = './DNABank-v3.csv'
    # start here
    if arg.query_file is None and arg.query_list is None:
        raise Exception('Wrong input!')
    query()
    # end


if __name__ == '__main__':
    main()
