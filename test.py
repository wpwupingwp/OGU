#!/usr/bin/python3

from functools import wraps
from timeit import default_timer as timer
from Bio import AlignIO
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
def test_np(a):
    length = len(a[0])
    for i in range(length):
        c = a[:, i]
        n = (c=='A').sum()
        n = (c=='T').sum()
        n = (c=='C').sum()
        n = (c=='G').sum()
        n = (c=='-').sum()


@print_time
def convert(a, order='C'):
    a = np.array([list(i) for i in a], order=order)
    return a


@print_time
def test_bio(a):
    length = len(a[0])
    for i in range(length):
        c = a[:, i]
        n = c.count('A')
        n = c.count('T')
        n = c.count('C')
        n = c.count('G')
        n = c.count('-')


@print_time
def read():
    a = AlignIO.read('./aln', 'fasta')
    return a


@print_time
def main():
    """The result shows that convert+numpy is abouth twice faster than
    Bio.AlignIO iterate. And the order 'F' is faster than 'C' order(~1.4
    times).
    """
    # start here
    a = read()
    test_bio(a)
    a = convert(a)
    test_np(a)
    a = read()
    a = convert(a, order='F')
    test_np(a)
    # end


if __name__ == '__main__':
    main()
