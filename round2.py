#!/usr/bin/python3

from __future__ import print_function
import glob
import pandas
from scipy.stats import linregress


def get_raw_data(filenames, data):
    """Here it uses pandas.read_excel to get contents from all xls files in
    present directory. Since pandas uses numpy.int64 to store data, it may
    cause a problem to manipulate huge numbers(larger than 2**63).  
    """
    for name in filenames:
        sheet = pandas.read_excel(name)
# drop '.xls' from name
        data[name[:5]] = sheet


def initiate_sample_data(raw, sample):
    """This function will get every sample's data and its two references 
    data.
    table:
        >>> a.index
        Index(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'], dtype='object')
        >>> a.columns
        Index([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 'Unnamed: 12'],
              dtype='object')
    sample name then will looks like:
        A33-11D
    """
    for library in ['A', 'B']:
        for plate in range(56):
            if library == 'B' and plate > 25:
                continue
            for idx in 'ABCDEFGH':
                for col in ['02', '03', '04', '05', '06',
                            '07', '08', '09', '10', '11']:
                    name = ''.join([
                        library,
                        '{:{fill}2d}'.format(plate+1, fill='0'),
                        '-',
                        col,
                        idx
                    ])
                    sample[name] = { 
                        'raw': [0, 0, 0, 0, 0, 0],
                        'ref_1': [0, 0, 0, 0, 0, 0],
                        'ref_2': [0, 0, 0, 0, 0, 0]
                    }


def get_sample_data(raw_data, sample):
    """It is better to check the original data to ensure this function works
    well.
    """
    for sheet_name, sheet in raw_data.items():
        time = int(sheet_name[-1])
        if time > 6:
            continue
        else:
            time = time - 1
        for idx in 'ABCDEFGH':
            for col in ['02', '03', '04', '05', '06', '07', '08', '09', '10', '11']:
                cell = ''.join([
                    sheet_name[0:-1],
                    col,
                    idx
                ])
                ref_1 = sheet[1][idx]
                ref_2 = sheet[12][idx]
                raw = sheet[int(col)][idx]
                sample[cell]['raw'][time] = raw
                sample[cell]['ref_1'][time] = ref_1
                sample[cell]['ref_2'][time] = ref_2


def analyse(sample_raw_data, analysis):
    """This function only use the first five points.
    """
    x = [0, 60, 120, 180, 240]
    for name, data in sample_raw_data.items():
        id = name
        item = [id, 0, 0]
        raw = data['raw'][:5]
        ref_1 = data['ref_1'][:5]
        ref_2 = data['ref_2'][:5]
        if 'OVRFLW' in raw:
            continue
        slope, intercept, r_value, _, _ = linregress(x, raw)
        r_square = r_value ** 2
        item.extend([
            slope,
            intercept,
            r_square
        ])
        slope, intercept, r_value, _, _ = linregress(x, ref_1)
        r_square = r_value ** 2
        item.extend([
            slope,
            intercept,
            r_square
        ])
        slope, intercept, r_value, _, _ = linregress(x, ref_2)
        r_square = r_value ** 2
        item.extend([
            slope,
            intercept,
            r_square
        ])
        item.extend(raw)
        item.extend(ref_1)
        item.extend(ref_2)
        item[1] = item[3] / item[6]
        item[2] = item[3] / item[9]
        analysis.append(item)


def output(analysis):
    """Output csv format.
    """
    with open('result.csv', 'w') as out:
        for line in analysis:
            line_out = [str(i) for i in line]
            out.write(','.join(line_out))
            out.write('\n')


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
    name_list = glob.glob('*-*')
    raw_data = dict()
    sample_raw_data = dict()
    sample = dict()
    analysis = list()
    analysis = [[
        'id', 
        'cell',
        'slope_of_slope',
        'fold_1', 'fold_2', 'fold_3',
        'slope_50uM', 'slope_20uM', 'slope_10uM', 
        'slope_ref1', 'slope_ref2', 'slope_ref3',
        'intercept_50uM', 'intercept_20uM', 'intercept_10uM',
        'intercept_ref1', 'intercept_ref2', 'intercept_ref3',
        'r^2_50uM', 'r^2_20uM', 'r^2_10uM',
        'r^2_ref1', 'r^2_ref2', 'r^2_ref3',
        'raw_1', 'raw_2', 'raw_3', 'raw_4', 'raw_5', 'raw_6', 'raw_7',
        'ref1_1', 'ref1_2', 'ref1_3', 'ref1_4', 'ref1_5', 'ref1_6', 'ref1_7', 
        'ref2_1', 'ref2_2', 'ref2_3', 'ref2_4', 'ref2_5', 'ref2_6', 'ref2_7', 
    ]]
    get_raw_data(name_list, raw_data)
    initiate_sample_data(raw_data, sample_raw_data)
    get_sample_data(raw_data, sample_raw_data)
    analyse(sample_raw_data, analysis)
    output(analysis)

if __name__ == '__main__':
    main()
