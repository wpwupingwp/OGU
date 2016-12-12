#!/usr/bin/python3

import argparse
from matplotlib import pyplot as plt
from matplotlib.markers import MarkerStyle
from random import choice


def get_data(data_file):
    """
    Hyposis that there is no SD bar for x axis.
    """
    unit = dict()
    y = dict()
    y_sd = dict()
    with open(data_file, 'r') as raw:
        for line in raw:
            line = line.strip()
            line = line.split(sep=arg.split)
            if line[0].startswith('x_unit'):
                unit['x'] = line[0].split(sep='=')[1]
                continue
            if line[0].startswith('y_unit'):
                unit['y'] = line[0].split(sep='=')[1]
                continue
            if line[0].startswith('x'):
                x = line[1:]
                x = [float(i) for i in x]
                continue
            if line[0].startswith('y'):
                y_id, label = line[0].split(sep=':')
                if label != 'sd':
                    y[y_id] = [float(i) for i in line[1:]]
                else:
                    y_sd[y_id] = [float(i) for i in line[1:]]
    return unit, x, y, y_sd


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
    unit, x, y, y_sd = data
    plt.xlabel(unit['x'])
    plt.ylabel(unit['y'])
    plt.plot(x, y['y2'])
    for i in y.keys():
        plt.errorbar(x, y[i], yerr=y_sd[i], fmt='k-'+choice(markers), label=i)
    plt.show()
    # end


if __name__ == '__main__':
    main()
