#!/usr/bin/python3

import argparse
from timeit import default_timer as timer


def main():
    """docstring
    """
    start_time = timer()
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('--path', default='./',
                        help='target path, default is "./"')
    parser.print_help()
    arg = parser.parse_args()
    # start here
    print(arg)
    # end
    end_time = timer()
    print('Cost {:.3f}s.\n'.format(end_time-start_time))

if __name__ == '__main__':
    main()
