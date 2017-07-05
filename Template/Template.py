#!/usr/bin/python3

import argparse
import os
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
def function():
    pass


def parse_args():
    arg = argparse.ArgumentParser(description=main.__doc__)
    arg.add_argument('-o', '--out', default='out',
                     help='output directory')
    arg.print_help()
    return arg.parse_args()


def main():
    """docstring
    """
    arg = parse_args()
    # start here
    if not os.path.exists(arg.out):
        os.mkdir(arg.out)
    function()
    # end


if __name__ == '__main__':
    main()
