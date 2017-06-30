#!/usr/bin/python3

from glob import glob
from multiprocessing import Pool, cpu_count
from subprocess import run
from timeit import default_timer as timer
import argparse


def function(parameters):
    string, filename = parameters
    string = string.replace('%i', '{}')
    string = string.format(filename)
    print(string)
    run(string, shell=True)
    return


def main():
    start = timer()
    arg = argparse.ArgumentParser()
    arg.add_argument('command',
                     help='command to run, use "$i" to stand for file name')
    arg.add_argument('files')
    arg.add_argument('-cpu', default=cpu_count()-1, help='CPU cores to use')
    arg = arg.parse_args()
    print(arg)

    files = glob(arg.files)
    tasks = [(arg.command, i) for i in files]
    pool = Pool(arg.cpu)
    result = pool.map(function, tasks)
    pool.terminate()
    pool.join()
    print(result)
    end = timer()

    print('\nFinished with {0:.3f}s.\n'.format(end-start))


if __name__ == '__main__':
    main()
