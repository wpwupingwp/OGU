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
#drop '.xls' from name
        data[name[:5]] = sheet

def get_sample_data(raw, sample):
    """This function will get every sample's data and its two references 
    data.
    """
#Initiate sample, aka sample_raw_data
    default = {
        'raw':[0,0,0,0,0,0],
        'ref_1':[0,0,0,0,0,0],
        'ref_2':[0,0,0,0,0,0]
    }
    for library in ['A', 'B']:
        for plate in range(56):
            if library == 'B' and plate>25:
                continue
            for line in 'ABCDEFGH':
                for row in range(12):
                    name = ''.join([
                        library,
                        '{:{fill}2d}'.format(plate+1,fill='0'),
                        '-',
                        line,
                        '{:{fill}2d}'.format(row+1,fill='0'),
                    ])
                    sample[name] = default
                #value = sheet[row][line]


    print(sample)
    raise ValueError
    for sheet_name in raw.keys():
        sheet = raw[sheet_name]
        for row in sheet.columns:
#row 1 and row 12 are references
            if row in [1, 12]:
                continue
            else:
                for line in sheet.index:
                    name = ''.join([
                        library,
                        plate,
                        '-',
                        sheet.index,
                        sheet.columns
                    ])
                    value = sheet[row][line]


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
    name_list = glob.glob('*.xls')
    get_raw_data(name_list, raw_data)
    sample_raw_data = dict()
    get_sample_data(raw_data, sample_raw_data)
    print(sample_raw_data)


if __name__ == '__main__':
    main()
