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


# @print_time
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
                record.append(line[:-1].upper())
        data.append([record[0], ''.join(record[1:])])
    data = data[1:]
    return data


# @print_time
def convert(old):
    # order 'F' is a bit faster than 'C'
    # name = np.array([[i[0]] for i in old], dtype=np.bytes_)
    # seq = np.array([list(i[1]) for i in old], dtype=np.bytes_)
    # new = np.hstack((name, seq))
    new = np.array([list(i[1]) for i in old], dtype=np.bytes_)
    rows, columns = new.shape
    return new, rows, columns


# @print_time
def count(alignment, rows, columns):
    # skip sequence id column
    data = [[rows, columns]]
    for index in range(1, columns):
        column = alignment[:, [index]]
        a = (column == b'A').sum()
        t = (column == b'T').sum()
        c = (column == b'C').sum()
        g = (column == b'G').sum()
        gap = rows - a - t - c - g
        data.append([a, t, c, g, gap])
    return data


# @print_time
def find_most(data, cutoff, gap_cutoff):
    rows, columns = data[0]
    data = data[1:]
    gap_cutoff = rows * gap_cutoff
    cutoff = rows * cutoff
    most = [['location', 'base', 'count']]
    for location, column in enumerate(data, 1):
        a, t, c, g, gap = column
        if gap >= gap_cutoff:
            continue
        if a >= cutoff:
            base = 'A'
            count = a
        elif t >= cutoff:
            base = 'T'
            count = t
        elif c >= cutoff:
            base = 'C'
            count = c
        elif g >= cutoff:
            base = 'G'
            count = g
        else:
            base = '?'
            count = gap
        most.append([location, base, count])
    return most[1:]


# @print_time
def find_continuous(most, window):
    continuous = list()
    fragment = list()
    most = [i for i in most if i[1] not in ('?', 'N')]
    for index, value in enumerate(most):
        fragment.append(value)
        location, *_ = value
        try:
            location_next, *_ = most[index+1]
        except:
            fragment.append(value)
            break
        step = location_next - location
        if step > 1:
            continuous.append(fragment)
            fragment = list()
    return continuous


# @print_time
def find_primer(continuous, most, window):
    for i in continuous:
        start = i[0][0]
        end = i[-1][0]
        seq = ''.join([j[1] for j in i])
        if end - start >= window:
            print(start, end, end-start, seq, sep='\t')


# @print_time
def parse_args():
    arg = argparse.ArgumentParser(description=main.__doc__)
    arg.add_argument('input', help='input alignment file')
    arg.add_argument('-c', '--cutoff', type=float, default=1.0,
                     help='minium percent to keep')
    arg.add_argument('-g', '--gap_cutoff', type=float, default=0.5,
                     help='maximum percent for gap to cutoff')
    arg.add_argument('-o', '--output', default='new.fasta')
    arg.add_argument('-w', '--window', type=int, default=20,
                     help='swip window width')
    # arg.print_help()
    return arg.parse_args()


# @print_time
def main():
    """Use ~4gb. Try to reduce.
    """
    arg = parse_args()
    alignment = read(arg.input)
    new, rows, columns = convert(alignment)
    print('{} has {} sequences with {} width.'.format(arg.input, rows,
                                                      columns))
    count_data = count(new, rows, columns)
    most = find_most(count_data, arg.cutoff, arg.gap_cutoff)
    continuous = find_continuous(most, arg.window)
    find_primer(continuous, most, arg.window)


if __name__ == '__main__':
    main()
