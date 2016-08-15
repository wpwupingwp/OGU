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
    length = len(raw)
    fields = len(raw[0])
    n = 1
    while n < length - 1:
        line_present = raw[n]
        line_after = raw[n+1]
        offset = 0
        merge = [set([i]) for i in line_present]
        while line_present[key] == line_after[key]:
            offset += 1
            for i, field in enumerate(line_after):
                merge[i].add(field)
            line_present = raw[n+offset]
            line_after = raw[n+1+offset]
        merge = ['/'.join(j) for j in merge]
        new.append(merge)
        n = n + offset + 1

    with open('output.csv', 'w', newline='', encoding='utf8') as output_file:
        writer = csv.writer(output_file)
        writer.writerows(new)

    print('Cost {:.3f}s.\n'.format(process_time()))

if __name__ == '__main__':
    main()
