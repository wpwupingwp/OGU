#!/usr/bin/python3

from time import process_time
import argparse
import csv

def main():
    '''This program will merge duplicate line in given table according to the key.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='file to be processed, must be csv')
    parser.add_argument('column', help='which column set as key, start at 0')
    arg = parser.parse_args()

    with open(arg.file, 'r', encoding='utf8') as f:
        reader = csv.reader(f)
        raw = [line for line in reader]
    key = int(arg.column)
    new = list()
    new.append(raw[0])

    for n,line_present in enumerate(raw):
        if n == 1:
            continue
        line_before = raw[n-1]
        if line_present[key] == line_before[key]:
            for m, field in enumerate(line_present):
                ref = line_before[m]
                if field != ref:
                    if field == '':
                        field = ref
                    else:
                        field = ref + '/' + field
            new.pop(-1)
        new.append(line_present)

    with open('output.csv', 'w', newline='', encoding='utf8') as output_file:
        writer = csv.writer(output_file)
        writer.writerows(new)

    print('Cost {:.3f}s.\n'.format(process_time()))
    done = input('Press any key to continue ...')

if __name__ == '__main__':
    main()
