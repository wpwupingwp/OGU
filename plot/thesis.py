#!/usr/bin/python3

import argparse
from matplotlib import pyplot as plt
from matplotlib.markers import MarkerStyle
from random import choice


def convert(line, target='float'):
    if target == 'float':
        return [float(i) for i in line]
    elif target == 'int':
        return [int(i) for i in line]


def get_data(data_file):
    """
    Hyposis that there is no SD bar for x axis.
    """
    axis_unit = dict()
    y = dict()
    with open(data_file, 'r') as raw:
        raw = raw.readlines()
        raw = [i.strip() for i in raw]
        raw = [i.split(sep=arg.split) for i in raw]
    for n, line in enumerate(raw):
        label = line[0]
        unit = line[1]
        value = convert(line[2:])
        if label.startswith('x_value'):
            axis_unit['x'] = unit
            x = convert(value)
            continue
        elif label.startswith('#'):
            continue
        else:
            axis_unit['y'] = unit
        last_line = raw[n-1]
        if label == last_line[0]:
            y[last_line[0]] = (convert(last_line[2:]), convert(value))
        else:
            y[label] = (convert(value), list())
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
    parser.add_argument('--path', default='./',
                        help='target path, default is "./"')
    parser.add_argument('data', default='a.txt')
    parser.add_argument('-t', '--type', dest='type',
                        choices=['dot', 'bar', 'line, dotline'],
                        default='line')
    parser.add_argument('-s', '--split', default=' ', type=str)
    parser.add_argument('-o', '--output', dest='output',
                        default='output')
    global arg
    arg = parser.parse_args()
    markers = MarkerStyle.filled_markers
    data = get_data(arg.data)
    unit, x, y = data
    plt.xlabel(unit['x'], fontsize=16)
    plt.ylabel(unit['y'], fontsize=16)
    if arg.type == 'line':
        fmt = 'k-'
    elif arg.type == 'dot':
        fmt = 'ko'
    elif arg.type == 'dot_line':
        fmt = 'ko-'
    for i in y.keys():
        if len(y[i][1]) == 0:
            # 'k' means black
            plt.plot(x, y[i][0], fmt, label=i)
        else:
            fmt = 'k-'
            plt.errorbar(x, y[i][0], yerr=y[i][1],
                         fmt=fmt+choice(markers), label=i)
    if len(y) > 1:
        plt.legend(loc='best')
    plt.savefig(arg.output+'.png')
    plt.show()


if __name__ == '__main__':
    main()
