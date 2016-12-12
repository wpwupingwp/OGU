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
        else:
            axis_unit['y'] = unit
        last_line = raw[n-1]
        if label == last_line[0]:
            y[last_line[0]] = (convert(value), convert(last_line[2:]))
        else:
            y[label] = (convert(value), list())
    return axis_unit, x, y


def main():
    """docstring
    """
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('--path', default='./',
                        help='target path, default is "./"')
    parser.add_argument('data', default='a.txt')
    parser.add_argument('-t', '--type', dest='type',
                        choices=['dot', 'bar', 'line'])
    parser.add_argument('-s', '--split', default=' ', type=str)
    global arg
    arg = parser.parse_args()
    # start here
    data = get_data(arg.data)
    markers = MarkerStyle.filled_markers
    unit, x, y = data
    print(data)
    plt.xlabel(unit['x'])
    plt.ylabel(unit['y'])
    for i in y.keys():
        if len(y[i][1]) == 0:
            plt.plot(x, y[i], fmt='k-'+choice(markers), label=i)
        else:
            plt.errorbar(x, y[i][0], yerr=y[i][1],
                         fmt='k-'+choice(markers), label=i)
    plt.legend(loc='best')
    plt.show()
    # end


if __name__ == '__main__':
    main()
