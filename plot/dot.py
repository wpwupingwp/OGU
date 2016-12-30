#!/usr/bin/python3

import argparse
from matplotlib import pyplot as plt
from scipy.stats.stats import linregress

import matplotlib
matplotlib.rcParams['lines.linewidth'] = 1.5
matplotlib.rcParams['axes.linewidth'] = 1.5
matplotlib.rcParams['axes.labelsize'] = 16
matplotlib.rcParams['font.size'] = 10


def convert(line, target='float'):
    if target == 'float':
        return [float(i) for i in line]
    elif target == 'int':
        return [int(i) for i in line]
    elif target == 'str':
        return [str(i) for i in line]


def get_data(data_file):
    """
    Hyposis that there is no SD bar for x axis.
    """
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
        value = convert(line[2:])
        if len(value) == 0:
            continue
        if label.startswith('x_value'):
            axis_unit['x'] = unit
            x = convert(value)
            continue
        elif label.startswith('#'):
            continue
        else:
            axis_unit['y'] = unit
        last_line = clean_data[n-1]
        if label == last_line[0]:
            y[last_line[0]] = [convert(last_line[2:]), convert(value)]
        else:
            y[label] = [convert(value), list()]
    return axis_unit, x, y


def draw_fit(x, y):
    range = arg.regression
    x = x[:range]
    y = y[:range]
    slope, intercept, r_value, *_ = linregress(x, y)
    text = r'$f(x)={0:.3f}x+{1:.3f}, R^2={2:.3f}$'.format(
        slope, intercept, r_value**2)
    fit = [slope*i+intercept for i in x]
    plt.plot(x[:arg.regression], fit, 'k--')
    plt.annotate(text, xy=(x[-1], y[-2]))


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
    parser.add_argument('-s', '--split', default=' ', type=str)
    parser.add_argument('-o', '--output', default='output')
    parser.add_argument('-r', '--regression', type=int, default=0)
    global arg
    arg = parser.parse_args()
    data = get_data(arg.data)
    unit, x, y = data
    plt.xlabel(unit['x'])
    plt.ylabel(unit['y'])
    markers = ('o', 'v', '^', '<', '>', 's', '8', 'p')
    # line format
    fmt = '-'
    matplotlib.rcParams['lines.linewidth'] = 2
    matplotlib.rcParams['axes.linewidth'] = 2
    for i, marker in zip(y.keys(), markers):
        if len(y[i][1]) == 0:
            y[i][1] = [0] * len(y[i][0])
        plt.errorbar(x, y[i][0], fmt=fmt+marker, markeredgecolor='none',
                     yerr=y[i][1], label=i)
        if arg.regression != 0:
            draw_fit(x, y[i][0])
    if len(y) > 1:
        plt.legend(loc='best')
    plt.savefig(arg.output+'.png')
    plt.show()


if __name__ == '__main__':
    main()
