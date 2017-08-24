#!/usr/bin/python3

from functools import wraps
from timeit import default_timer as timer
import argparse
import numpy as np


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
def read(fasta):
    data = list()
    record = ['id', 'seq']
    with open(fasta, 'r') as raw:
        for line in raw:
            if line.startswith('>'):
                data.append([record[0], ''.join(record[1:])])
                name = line[1:-1]
                record = [name, ]
            else:
                record.append(line[:-1])
        data.append([record[0], ''.join(record[1:])])
    data = data[1:]
    return data


@print_time
def convert(old):
    # order 'F' is a bit faster than 'C'
    name = np.array([[i[0]] for i in old], dtype=np.bytes_)
    seq = np.array([list(i[1]) for i in old], dtype=np.bytes_)
    new = np.hstack((name, seq))
    return new, new.shape


@print_time
def count(alignment, columns, rows):
    # skip sequence id column
    data = [[columns, rows]]
    for index in range(1, rows):
        column = alignment[:, [index]]
        a = (column == b'A').sum()
        t = (column == b'T').sum()
        c = (column == b'C').sum()
        g = (column == b'G').sum()
        gap = columns - a - t - c - g
        data.append([a, t, c, g, gap])
    return data


def find_conservative(data, arg):
    rows, columns = data[0]
    data = data[1:]
    gap_cutoff =  rows * arg.gap_cutoff
    cutoff = rows * arg.cutoff
    most = [['location', 'base', 'count']]
    for index, column in enumerate(data, 1):
        a, t, c, g, gap = column
        if gap >= gap_cutoff:
            continue
        if a >= cutoff:
            base = 'A'
            count = a
        elif a >= cutoff:
            base = 'T'
            count = t
        elif a >= cutoff:
            base = 'C'
            count = c
        elif a >= cutoff:
            base = 'G'
            count = g
        else:
            base = '?'
            count = gap
        most.append([index, base, count])
    for i in most:
        print(*i)
    print(data)




@print_time
def write(data, output_file):
    with open(output_file, 'wb') as output:
        for i in data:
            seq_id, seq = i[0], i[1:]
            seq = b''.join(seq)
            output.write(b'>'+seq_id+b'\n'+seq+b'\n')


@print_time
def parse_args():
    arg = argparse.ArgumentParser(description=main.__doc__)
    arg.add_argument('input', help='input alignment file')
    arg.add_argument('-c', '--cutoff', type=float, default=0.99,
                     help='minium percent to keep')
    arg.add_argument('-g', '--gap_cutoff', type=float, default=0.5,
                     help='maximum percent for gap to cutoff')
    arg.add_argument('-o', '--output', default='new.fasta')
    arg.add_argument('-w', '--window', type=int, default=20,
                     help='swip window width')
    arg.print_help()
    return arg.parse_args()


@print_time
def main():
    """Use ~4gb. Try to reduce.
    """
    arg = parse_args()
    alignment = read(arg.input)
    new, shape = convert(alignment)
    count_data = count(new, *shape)
    find_conservative(count_data, arg)


if __name__ == '__main__':
    main()
