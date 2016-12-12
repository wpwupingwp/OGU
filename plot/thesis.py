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
def function():
    print('ok')
    pass


def main():
    """docstring
    """
    start_time = timer()
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('--path', default='./',
                        help='target path, default is "./"')
    parser.add_argument('data', default='a.txt')
    parser.add_argument('-t','--type', dest='type', choices=['dot', 'bar', 'line'])
    parser.print_help()
    arg = parser.parse_args()
    # start here
    print(vars(arg))
    function()
    # end
    end_time = timer()
    print('Cost {:.3f}s.\n'.format(end_time-start_time))

if __name__ == '__main__':
    main()
