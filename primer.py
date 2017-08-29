#!/usr/bin/python3

from functools import wraps
from collections import defaultdict
from timeit import default_timer as timer
import argparse
import numpy as np
from Bio.Data. IUPACData import ambiguous_dna_values


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


def get_ambiguous_dict():
    data = ambiguous_dna_values
    data = dict(zip(data.values(), data.keys()))
    # 2:{'AC': ['M',}
    data_with_len = defaultdict(lambda: dict())
    for key in data:
        data_with_len[len(key)][key] = data[key]
    return data_with_len


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
                record.append(line[:-1].upper())
        data.append([record[0], ''.join(record[1:])])
    data = data[1:]
    return data


@print_time
def convert(old):
    # order 'F' is a bit faster than 'C'
    # name = np.array([[i[0]] for i in old], dtype=np.bytes_)
    # seq = np.array([list(i[1]) for i in old], dtype=np.bytes_)
    # new = np.hstack((name, seq))
    new = np.array([list(i[1]) for i in old], dtype=np.bytes_, order='F')
    rows, columns = new.shape
    return new, rows, columns


@print_time
def count(alignment, rows, columns):
    # skip sequence id column
    data = [[rows, columns]]
    for index in range(columns):
        column = alignment[:, [index]]
        a = (column == b'A').sum()
        t = (column == b'T').sum()
        c = (column == b'C').sum()
        g = (column == b'G').sum()
        n = (column == b'N').sum()
        question = (column == b'?').sum()
        gap = (column == b'-').sum()
        # is it necessary to count 'N' '-' and '?' ?
        other = rows - a - t - c - g - n - question - gap
        data.append([a, t, c, g, n, question, gap, other])
    # make sure rows and columns does not mixed
    assert len(data) == columns + 1
    return data


@print_time
def find_most(data, cutoff, gap_cutoff):
    most = [['location', 'base', 'count']]
    rows, columns = data[0]
    data = data[1:]
    gap_cutoff = rows * gap_cutoff
    cutoff = rows * cutoff
    ambiguous_dict = get_ambiguous_dict()

    def run():
        for location, column in enumerate(data, 1):
            finish = False
            value = dict(zip(list('ATCGN?-X'), column))
            base = 'N'

            sum_gap = sum([value['?'], value['-'], value['X']])
            if sum_gap >= gap_cutoff:
                base = '-'
                count = sum_gap
                yield [location, base, count]
                continue
            # 1 2 3 4
            for length in ambiguous_dict:
                if finish:
                    break
                for key in ambiguous_dict[length]:
                    if finish:
                        break
                    count = 0
                    for letter in list(key):
                        if finish:
                            break
                        count += value[letter]
                        if count >= cutoff:
                            base = ambiguous_dict[length][key]
                            finish = True
                            yield [location, base, count]
    for i in run():
        most.append(i)
    return most[1:]


def print_consensus(data):
    i = [i[0] for i in data]
    seq = [i[1] for i in data]
    num = [i[2] for i in data]
    print('{:>5} {}> {:<5}'.format(i[0], '-'*5*(len(i)-2), i[-1]))
    for _ in seq:
        print('{:>5}'.format(_), end='')
    print()
    for _ in seq:
        print('{:>5}'.format('|'), end='')
    print()
    for _ in num:
        print('{:>5}'.format(_), end='')
    print()



@print_time
def find_continuous(most, window):
    continuous = list()
    fragment = list()
    most = [i for i in most if i[1] not in ('N', '-')]
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


@print_time
def find_primer(continuous, most, window):
    for i in continuous:
        if len(i) >= window:
            print_consensus(i)


@print_time
def parse_args():
    arg = argparse.ArgumentParser(description=main.__doc__)
    arg.add_argument('input', help='input alignment file')
    arg.add_argument('-c', '--cutoff', type=float, default=1.0,
                     help='minium percent to keep')
    arg.add_argument('-g', '--gap_cutoff', type=float, default=0.5,
                     help='maximum percent for gap to cutoff')
    arg.add_argument('-o', '--output', default='new.fasta')
    arg.add_argument('-w', '--window', type=int, default=18,
                     help='swip window width')
    # arg.print_help()
    return arg.parse_args()


@print_time
def main():
    arg = parse_args()
    raw_alignment = read(arg.input)
    new, rows, columns = convert(raw_alignment)
    count_data = count(new, rows, columns)
    most = find_most(count_data, arg.cutoff, arg.gap_cutoff)
    continuous = find_continuous(most, arg.window)
    find_primer(continuous, most, arg.window)


if __name__ == '__main__':
    main()
