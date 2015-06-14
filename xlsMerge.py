#!/usr/bin/python3

import glob
import pandas

def get_data(filenames, data):
    """Here it uses pandas.read_excel to get contents from all xls files in
    present directory. Since pandas uses numpy.int64 to store data, it may
    cause a problem to manipulate huge numbers(larger than 2**63).  
    """
    for name in filenames:
        sheet = pandas.read_excel(name)
        data[name] = sheet

def analyse(raw, edited):
    """Library: A B
       Plate: A 56, B 26
       Times: 6 
            There are A28 and others which have 7 or 10 tables. In order to
            simplify the program, here it just omit tables from 7 to 10.
    """
    for name in raw.keys():
        print(name[:5],type(raw[name]))
        library = name[0]
        plate = name[1:3]
        times = name[4]


def main():
    """It uses glob.glob to get names of all xls files. Hence it should be 
    run in the directory which contains all xls files.
    All xls filename should follow this format:
     A00-0.xls
     01234
    """
    raw_data = dict()
    analysed = list()
    name_list = glob.glob('*.xls')
    get_data(name_list, raw_data)
    analyse(raw_data, analysed)


if __name__ == '__main__':
    main()
