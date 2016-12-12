#!/usr/bin/python3

import argparse
import seaborn as s


def set_seaborn():
    pass


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
                x_label = line[0].split(sep=':')[0]
                x = line[1:]
                continue
            if line[0].startswith('y'):
                y_id, label = line[0].split(sep=':')
                if label != 'sd':
                    y[y_id] = line[1:]
                else:
                    y_sd = line[1:]
    return unit, x, x_label, y, y_sd


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
    set_seaborn()
    data = get_data(arg.data)
    for i in data:
        print(i)
    unit, x, x_label, y, y_sd = data
    # end


if __name__ == '__main__':
    main()
