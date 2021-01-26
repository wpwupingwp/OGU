#!/usr/bin/python3

from sys import argv
from BarcodeFinder.bf import bf_main
from BarcodeFinder.utils import get_all_third_party

def main():
    if argv[-1] == 'init':
        get_all_third_party()
    else:
        bf_main()


if __name__ == '__main__':
    main()