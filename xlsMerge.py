#!/usr/bin/python3

import glob
import pandas

def get_raw_data(filenames, data):
    """Here it uses pandas.read_excel to get contents from all xls files in
    present directory. Since pandas uses numpy.int64 to store data, it may
    cause a problem to manipulate huge numbers(larger than 2**63).  
    """
    for name in filenames:
        sheet = pandas.read_excel(name)
        data[name] = sheet

def get_sample_data(raw, sample):
    """This function will get every sample's data and its two references 
    data.
    """
    for library in ['A', 'B']:
        for plate in range(56):
            for times in range(6):
                sheet_name = ''.join([
                    library,
                    '{:{fill}2d}'.format(plate+1,fill='0'),
                    '-',
                    str(times+1)
                ])
                sheet = raw[sheet_name]
                print(sheet.index,shet.columns)




def main():
    """It uses glob.glob to get names of all xls files. Hence it should be 
    run in the directory which contains all xls files.
    All xls filename should follow this format:
     A00-0.xls
     01234

    Library: A B
       Plate: A 56, B 26
       Times: 6 
            There are A28 and others which have 7 or 10 tables. In order to
            simplify the program, here it just omit tables from 7 to 10.
    """
    raw_data = dict()
    sample_raw_data = {
        'raw':[1, 2, 3, 4, 5, 6],
        'ref_1':[1, 2, 3, 4, 5, 6],
        'ref_2':[1, 2, 3, 4, 5, 6]
    }

    name_list = glob.glob('*.xls')
    get_raw_data(name_list, raw_data)
    get_sample_data(raw_data, sample_raw_data)


if __name__ == '__main__':
    main()
