#!/usr/bin/python3

import glob
import pandas

def get_data(filenames, data):
    for name in filenames:
        sheet = pandas.read_excel(name)
        data.append([name, sheet])

def main():
    raw_data = list()
    name_list = glob.glob('*.xls')
    get_data(name_list, raw_data)
    print(raw_data)


if __name__ == '__main__':
    main()
