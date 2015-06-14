#!/usr/bin/python3

import glob
import pandas

def get_data(filenames, data):
    """Here it uses pandas.read_excel to get contents from all xls files in present directory. Since pandas uses numpy.int64 to store data, it may be a problem to manipulate huge numbers(larger than 2**63).
    """
    for name in filenames:
        sheet = pandas.read_excel(name)
        data.append([name, sheet])

def analyse():


def main():
    """It uses glob.glob to get names of all xls files. Hence it should be run in the directory which contains all xls files.
    """
    raw_data = list()
    name_list = glob.glob('*.xls')
    get_data(name_list, raw_data)
    analyse(raw_data)
    print(raw_data)


if __name__ == '__main__':
    main()
