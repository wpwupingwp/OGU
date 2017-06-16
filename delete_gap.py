#!/usr/bin/python3

import argparse
from functools import wraps
from timeit import default_timer as timer
from tempfile import mkdtemp
from os import path, mkdir


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
    tmp = mkdtemp()
    print(tmp)
    pass


def main():
    """docstring
    """
    parameters = argparse.ArgumentParser(description=main.__doc__)
    parameters.add_argument('--path', default='./',
                            help='target path, default is "./"')
    parameters.add_argument('-o', '--output', default='out',
                            help='output directory')
    parameters.print_help()
    arg = parameters.parse_args()
    # start here
    global tmp
    tmp = mkdtemp()
    print(vars(arg))
    if not path.exists(arg.output):
        mkdir(arg.output)
    function()
    # end


if __name__ == '__main__':
    main()
