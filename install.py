#!/usr/bin/python3

from json import load
from os import environ, pathsep, sep
from os.path import abspath
from platform import system
from shutil import unpack_archive, ReadError
from subprocess import run
from urllib import request


with open('url.json', 'r') as _:
    urls = load(_)
sys = system()
for software in urls[sys]:
    url = urls[software]
    filename = url.split('/')[-1]
    down = request.urlopen(url)
    if down.status == 200:
        with open(filename, 'wb') as out:
            out.write(down.read())
    else:
        print('Cannot download. Retry... (Ctrl+C to abort)')
        continue
    try:
        unpack_archive(filename)
    except ReadError:
        pass
if sys == 'Windows':
    environ['PATH'] = pathsep.join([abspath('mafft-win'),
                                    abspath('iqtree-1.6.7-Windows'+sep+'bin'),
                                    environ['PATH']])
    run('ncbi-blast-2.7.1+-win64.exe', shell=True)
elif sys == 'Linux':
    environ['PATH'] = pathsep.join([abspath('mafft-linux64'),
                                    abspath('iqtree-1.6.7-Linux'+sep+'bin'),
                                    abspath('ncbi-blast-2.7.1+'+sep+'bin'),
                                    environ['PATH']])
elif sys == 'Apple':
    pass
print(environ['PATH'])
