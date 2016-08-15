#!/usr/bin/python3

from time import process_time
import argparse
import csv

def main():
    parser = argparse.ArgumentParser(
        description='''This program will merge duplicate record in table.''')
    parser.add_argument('file', help='file to be processed, must be csv')
    parser.add_argument('column', help='which column set as key, start at 0')
    arg = parser.parse_args()

    with open(arg.file, 'r', encoding='utf8') as f:
        reader = csv.reader(f)
        raw = [line for line in reader]

    print('Cost {:.3f}s.\n'.format(process_time()))
    done = input('Press any key to continue ...')

if __name__ == '__main__':
    main()
