#!/usr/bin/python3

import argparse
from Bio import SeqIO
from glob import glob


def main():
    """This program will try to find out barcode to devide different species
    while ignore distinction among subspecies level."""
    parser = argparse.ArgumentParser( description=main.__doc__)
    parser.add_argument('--path', default='.', 
                        help='target path, default is present directory')
    parser.
    arg = parser.parse_args() 


if __name__ == '__main__':
    main()