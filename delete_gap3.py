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
    record = ['ID', 'Sequence']
    with open(fasta, 'r') as raw:
        for line in raw:
            if line.startswith('>'):
                record = [record[0], ''.join(record[1:])]
                data.append(record)
                record = [line[1:-1], ]
            else:
                record.append(line[:-1])
    data = data[1:]
    print(data[0])
    print(data[-1])
    return data


@print_time
def convert(old):
    # order 'F' is a bit faster than 'C'
    shape = (len(old), len(old[0]))
    print(shape)
    new = np.array([list(i) for i in old], order='F').reshape(shape)
    return new, new.shape


@print_time
def remove_gap(alignment, length, width):
    # get alignment head
    short = alignment[:, 0]
    for index in range(1, width):
        column = alignment[:, index]
        a = (column == 'A').sum()
        t = (column == 'T').sum()
        c = (column == 'C').sum()
        g = (column == 'G').sum()
        gap = length - a - t - c - g
        if gap == length:
            print('Empty in column {}'.format(index))
            print(gap, a, t, c, g, length)
            continue
        elif (a+1 == (length-1) or t+1 == (length-1) or c+1 == (length-1) or
              g+1 == (length-1)):
            print('Only one in column {}'.format(index))
            continue
        elif a == length or t == length or c == length or g == length:
            print('All same in column {}'.format(index))
            continue
        else:
            short = np.hstack((short, column))
            print(short.shape)
            continue
    return short, short.shape


@print_time
def write(data, output_file):
    with open(output_file, 'w') as output:
        for i in data:
            id, *seq = i
            output.write('>{}\n{}\n'.format(id, ''.join(seq)))


@print_time
def parse_args():
    arg = argparse.ArgumentParser(description=main.__doc__)
    arg.add_argument('input', help='input alignment file')
    arg.add_argument('-o', '--output', default='new.fasta')
    arg.print_help()
    return arg.parse_args()


@print_time
def main():
    """docstring
    """
    arg = parse_args()
    alignment = read(arg.input)
    new_alignment, shape = convert(alignment)
    after_delete, new_shape = remove_gap(new_alignment, *shape)
    write(after_delete, arg.output)
    print(new_shape, shape)
    print('Remove {} columns.'.format(new_shape[1]-shape[1]))


if __name__ == '__main__':
    main()
