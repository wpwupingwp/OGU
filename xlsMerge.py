#!/usr/bin/python3

import glob
import pandas

def get_data(filenames, data):
    """Here it uses pandas.read_excel to get contents from all xls files in
    present directory. Since pandas uses numpy.int64 to store data, it may be
    a problem to manipulate huge numbers(larger than 2**63).  
    """
    for name in filenames:
        sheet = pandas.read_excel(name)
        data.append([name, sheet])

def analyse(raw, edited):
    """Library: A B
       Plate: A 56, B 26
       Times: 6 7 9
       !rm A28-10.xls
    """
    for name, sheet in raw:
        print(name[:5],type(sheet))

def main():
    """It uses glob.glob to get names of all xls files. Hence it should be run
    in the directory which contains all xls files.
    All xls filename should follow this principle:
    A00-0.xls
    """
    raw_data = list()
    analysed = list()
    name_list = glob.glob('*.xls')
    get_data(name_list, raw_data)
    analyse(raw_data, analysed)


if __name__ == '__main__':
    main()
