#!/usr/bin/python3

from __future__ import print_function
import glob
import pandas
from scipy.stats import linregress


def get_raw_data(filenames, data):
    """It can accept xls format.
    """
    for name in filenames:
        sheet = pandas.read_excel(name)
        data[name] = sheet


def initiate_sample_data(sample):
    """This function will get every sample's data and its two references 
    data.
    table:
        >>> a.index
        Index(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'], dtype='object')
        >>> a.columns
        Index([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 'Unnamed: 12'],
              dtype='object')
    sample name then will looks like:
        P1-11D
    """
    for plate in '1234567':
        for idx in 'ABCDEFGH':
            for col in ['01', '12']:
                name = ''.join([
                    'P',
                    plate,
                    '-',
                    col,
                    idx
                ])
                sample[name] = { 
                    '50uM': [0, 0, 0, 0, 0, 0, 0],
                    '20uM': [0, 0, 0, 0, 0, 0, 0],
                    '10uM': [0, 0, 0, 0, 0, 0, 0],
                    'ref1': [0, 0, 0, 0, 0, 0, 0],
                    'ref2': [0, 0, 0, 0, 0, 0, 0],
                    'ref3': [0, 0, 0, 0, 0, 0, 0],
                    'ref4': [0, 0, 0, 0, 0, 0, 0],
                    'ref5': [0, 0, 0, 0, 0, 0, 0],
                    'ref6': [0, 0, 0, 0, 0, 0, 0]
                }


def get_sample_data(raw_data, sample):
    """In round 2, data in left plate and right plate are in different order:
            ref1 50uM 20uM 10uM ref2 ref3 50uM 20uM 10uM ref1 ref2 ref3
    But in Plate 1 and 2, there are some data does not follow this rule.
    """
    for sheet_name, sheet in raw_data.items():
        time = int(sheet_name[-1])-1
        plate = int(sheet_name[1])
        for idx in 'ABCDEFGH':
            for col in ['01', '12']:
                cell = ''.join([
                    sheet_name[0:-1],
                    col,
                    idx
                ])
                ref1 = sheet[1][idx]
                ref2 = sheet[12][idx]
                ref3 = ref1
                ref4 = ref2
                ref5 = ref1
                ref6 = ref2
                if plate >= 3:
                    ref3 = sheet[6][idx]
                    ref4 = sheet[11][idx]
                if plate >= 5:
                    ref5 = sheet[5][idx]
                    ref6 = sheet[10][idx]
                if plate == 7:
                    ref2 = ref1
                    ref4 = ref3
                    ref6 = ref5
                if col == '01':
                    fifty = sheet[2][idx]
                    twenty = sheet[3][idx]
                    ten = sheet[4][idx]
                else:
                    fifty = sheet[7][idx]
                    twenty = sheet[8][idx]
                    ten = sheet[9][idx]
                sample[cell]['50uM'][time] = fifty
                sample[cell]['20uM'][time] = twenty
                sample[cell]['10uM'][time] = ten
                sample[cell]['ref1'][time] = ref1
                sample[cell]['ref2'][time] = ref2
                sample[cell]['ref3'][time] = ref3
                sample[cell]['ref4'][time] = ref4
                sample[cell]['ref5'][time] = ref5
                sample[cell]['ref6'][time] = ref6


def analyse(sample_raw_data, analysis, id_list):
    """This function use all seven points.
    """
    x = [0, 60, 120, 180, 240, 300, 360]
    for name, data in sample_raw_data.items():
        item = [id_list[name],
                name, 
                0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0
                ]
        fifty = data['50uM']
        twenty = data['20uM']
        ten = data['10uM']
        if 'OVRFLW' in fifty:
            continue
        ref1 = data['ref1']
        ref2 = data['ref2']
        ref3 = data['ref3']
        ref4 = data['ref4']
        ref5 = data['ref5']
        ref6 = data['ref6']
        ref = list()
        for i in range(7):
            to_mean = [ref1[i], ref2[i], ref3[i], ref4[i], ref5[i], ref6[i]]
            ref.append(sum(to_mean)/6)

        slope, intercept, r_value, _, _ = linregress(x, fifty)
        item[5] = slope
        item[9] = intercept
        item[13] = r_value ** 2
        slope, intercept, r_value, _, _ = linregress(x, twenty)
        item[6] = slope
        item[10] = intercept
        item[14] = r_value ** 2
        slope, intercept, r_value, _, _ = linregress(x, ten)
        item[7] = slope
        item[11] = intercept
        item[15] = r_value ** 2
        slope, intercept, r_value, _, _ = linregress(x, ref)
        item[8] = slope
        item[12] = intercept
        item[16] = r_value ** 2

        item[2] = item[5] / item[8]
        item[3] = item[6] / item[8]
        item[4] = item[7] / item[8]

        item.extend(fifty)
        item.extend(twenty)
        item.extend(ten)
        item.extend(ref)
        item.extend(ref1)
        item.extend(ref2)
        item.extend(ref3)
        item.extend(ref4)
        item.extend(ref5)
        item.extend(ref6)
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
     P0-0.xls
    """
    id_list = open('round2.csv', 'r').read().split(sep='\n')
    id_list = [i.split()[::-1] for i in id_list]
    id_list.pop()
    id_list = dict(id_list)
    name_list = glob.glob('*-*')
    raw_data = dict()
    sample_raw_data = dict()
    analysis = [[
        'id', 
        'cell',
        'fold_50uM', 'fold_20uM', 'fold_10uM',
        'slope_50uM', 'slope_20uM', 'slope_10uM', 'slope_ref',
        'intercept_50uM', 'intercept_20uM', 'intercept_10uM', 'intercept_ref',
        'r^2_50uM', 'r^2_20uM', 'r^2_10uM', 'r^2_ref',
        '50uM_1', '50uM_2', '50uM_3', '50uM_4', '50uM_5', '50uM_6', '50uM_7',
        '20uM_1', '20uM_2', '20uM_3', '20uM_4', '20uM_5', '20uM_6', '20uM_7',
        '10uM_1', '10uM_2', '10uM_3', '10uM_4', '10uM_5', '10uM_6', '10uM_7',
        'ref_1', 'ref_2', 'ref_3', 'ref_4', 'ref_5', 'ref_6', 'ref_7', 
        'ref1_1', 'ref1_2', 'ref1_3', 'ref1_4', 'ref1_5', 'ref1_6', 'ref1_7', 
        'ref2_1', 'ref2_2', 'ref2_3', 'ref2_4', 'ref2_5', 'ref2_6', 'ref2_7', 
        'ref3_1', 'ref3_2', 'ref3_3', 'ref3_4', 'ref3_5', 'ref3_6', 'ref3_7', 
        'ref4_1', 'ref4_2', 'ref4_3', 'ref4_4', 'ref4_5', 'ref4_6', 'ref4_7', 
        'ref5_1', 'ref5_2', 'ref5_3', 'ref5_4', 'ref5_5', 'ref5_6', 'ref5_7', 
        'ref6_1', 'ref6_2', 'ref6_3', 'ref6_4', 'ref6_5', 'ref6_6', 'ref6_7', 
    ]]
    get_raw_data(name_list, raw_data)
    initiate_sample_data(sample_raw_data)
    get_sample_data(raw_data, sample_raw_data)
    analyse(sample_raw_data, analysis, id_list)
    output(analysis)

if __name__ == '__main__':
    main()
