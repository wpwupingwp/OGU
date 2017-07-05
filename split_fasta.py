#!/usr/bin/python3

import argparse
import os
from timeit import default_timer as timer


def parse_args():
    arg = argparse.ArgumentParser(description=main.__doc__)
    arg.add_argument('-i', '--input', required=True, help='input file')
    arg.add_argument('-s', '--size', type=int, default=100000,
                     help='how many sequences one file have')
    arg.add_argument('-o', '--out', default='out',
                     help='output directory')
    arg.print_help()
    return arg.parse_args()


def main():
    """docstring
    """
    start = timer()
    arg = parse_args()
    # start here
    if not os.path.exists(arg.out):
        os.mkdir(arg.out)
    function()
    # end
    end = timer()
    print('Cost {1:3f}s.\n'.format(end-start))


if __name__ == '__main__':
    main()
