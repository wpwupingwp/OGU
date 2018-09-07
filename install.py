#!/usr/bin/python3

from json import load
from os import environ, pathsep, sep
from os.path import abspath
from platform import system
from shutil import unpack_archive, ReadError
from subprocess import run
from urllib import request


def download(sys):
    with open('url.json', 'r') as _:
        urls = load(_)
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
    return True


sys = system()
if sys == 'Windows':
    download(sys)
    run('ncbi-blast-2.7.1+-win64.exe', shell=True)
    environ['PATH'] = pathsep.join([abspath('mafft-win'),
                                    abspath('iqtree-1.6.7-Windows'+sep+'bin'),
                                    environ['PATH']])
elif sys == 'Linux':
    ok = False
    for pack_mgr in ('apt', 'dnf', 'yum', 'pkg'):
        r = run('{} install ncbi-blast+ iqtree mafft'.format(pack_mgr),
                shell=True)
        if r.returncode == 0:
            ok = True
            break
    if not ok:
        download(sys)
        environ['PATH'] = pathsep.join([
            abspath('mafft-linux64'), abspath('iqtree-1.6.7-Linux'+sep+'bin'),
            abspath('ncbi-blast-2.7.1+'+sep+'bin'), environ['PATH']])
elif sys == 'Apple':
    download(sys)
print(environ['PATH'])
