#!/usr/bin/python3

import argparse
from random import choice
from matplotlib import pyplot as plt

import matplotlib
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.linewidth'] = 2
matplotlib.rcParams['axes.labelsize'] = 16
matplotlib.rcParams['font.size'] = 12


def convert(line, target='float'):
    if target == 'float':
        return [float(i) for i in line]
    elif target == 'int':
        return [int(i) for i in line]
    elif target == 'str':
        return [str(i) for i in line]


def get_data(data_file):
    clean_data = list()
    with open(data_file, 'r') as raw:
        for line in raw.readlines():
            if line.startswith('#'):
                continue
            line = line.strip()
            line = line.split(sep=arg.split)
            if len(line) < 2:
                continue
            clean_data.append(line)

    axis_unit = dict()
    y = dict()
    for n, line in enumerate(clean_data):
        label = line[0]
        unit = line[1]
        if len(line[2:]) == 0:
            continue
        if label.startswith('x_value'):
            axis_unit['x'] = unit
            x = convert(line[2:], 'str')
            continue
        elif label.startswith('#'):
            continue
        else:
            axis_unit['y'] = unit
        value = convert(line[2:])
        last_line = clean_data[n-1]
        if label == last_line[0]:
            y[last_line[0]] = [convert(last_line[2:]), convert(value)]
        else:
            y[label] = [convert(value), list()]
    return axis_unit, x, y


def main():
    """
    The input file should look like this:
    x_value unit 1 2 3
    ylabel unit 23 3 5
    y_sd unit 3 3 3
    The line starts with # will be ignored.
    """
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('data', default='a.txt')
    parser.add_argument('-s', '--split', default=' ', type=str)
    parser.add_argument('-o', '--output', default='output')
    parser.add_argument('-w', '--width', default='0.35', type=float)
    global arg
    arg = parser.parse_args()

    data = get_data(arg.data)
    unit, x, y = data
    fig, ax = plt.subplots()
    plt.xlabel(unit['x'])
    plt.ylabel(unit['y'])
    # format
    patterns = ('//', '\\', '', '/', '')
    index = list(range(len(x)))
    index_with_offset = [i+arg.width*(len(y)-1)/2 for i in index]
    plt.xticks(index_with_offset, x)
    for n, i in enumerate(y.keys()):
        index_i = [i+arg.width*n for i in index]
        if len(y[i][1]) == 0:
            y[i][1] = [0] * len(y[i][0])
        plt.bar(index_i, y[i][0], arg.width, color='white',
                yerr=y[i][1], ecolor='black', edgecolor='black',
                label=i, align='center', hatch=choice(patterns))
    if len(y) > 1:
        plt.legend(loc='best')
    fig.autofmt_xdate()
    plt.savefig(arg.output+'.png')
    plt.show()


if __name__ == '__main__':
    main()
