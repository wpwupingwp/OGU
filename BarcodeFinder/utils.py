#!/usr/bin/python3

import re
import logging

from platform import system
from os.path import abspath, basename, exists, splitext, pathsep, sep
from os.path import join as join_path
from os import (cpu_count, devnull, environ, mkdir, pathsep, remove,
                rename,)
from urllib.error import HTTPError
from urllib.request import urlopen
from shutil import unpack_archive, ReadError
from subprocess import run

# define logger
FMT = '%(asctime)s %(levelname)-8s %(message)s'
DATEFMT = '%H:%M:%S'
logging.basicConfig(format=FMT, datefmt=DATEFMT, level=logging.INFO)
log = logging.getLogger('barcodefinder')
try:
    import coloredlogs
    coloredlogs.install(level=logging.INFO, fmt=FMT, datefmt=DATEFMT)
except ImportError:
    pass


def plastid_rename():
    """
    Use name database.
    """
    pass


def safe_average(x):
    """
    Safe average.
    """
    if len(x) == 0:
        return 0
    else:
        return sum(x) / len(x)


def safe_path(old):
    """
    Remove illegal character in file path or name.
    """
    return re.sub(r'\W', '_', old)


def check_tools():
    """
    Check dependent software, if not found, try to install.
    Return original PATH.
    Return None if failed.
    """
    if exists('PATH.txt'):
        log.info('Reading configuration of previous run from PATH.txt.')
        with open('PATH.txt', 'r', encoding='utf-8') as path_file:
            exists_path = path_file.read().strip()
            environ['PATH'] = pathsep.join([environ['PATH'], exists_path])
    f = open(devnull, 'w', encoding='utf-8')
    installed = []
    # blast use different option style, have to use dict
    tools_cmd = {'MAFFT': 'mafft --version',
                 'IQTREE': 'iqtree --version',
                 'BLAST': 'makeblastdb -version'}
    for tools in tools_cmd:
        check = run(tools_cmd[tools], shell=True, stdout=f, stderr=f)
        # mafft --help return 0 or 1 in different version, use --version
        # instead
        if check.returncode != 0:
            log.warning('Cannot find {}.'.format(tools))
            install_path = deploy(tools)
            if install_path is None:
                log.error('Failed to install {}. Please try to manually '
                          'install it (See README.md).'.format(tools))
                return None
            installed.append(install_path)
    # do not edit original PATH
    to_add = pathsep.join(installed)
    original = str(environ['PATH'])
    environ['PATH'] = pathsep.join([original, to_add])
    if len(installed) != 0:
        log.info('Installation info of dependent software was written into '
                 'PATH.txt')
        with open('PATH.txt', 'w', encoding='utf-8') as path_out:
            path_out.write(to_add + '\n')
    else:
        log.info('All dependent software found.')
    f.close()
    return original


def download_software(url):
    """
    Download, return False if failed.
    http_proxy may affect this function.
    """
    filename = url.split('/')[-1]
    try:
        log.info('Downloading {}...'.format(filename))
        down = urlopen(url)
    except HTTPError:
        log.warning('Cannot download {}.'.format(filename))
        return False
    with open(filename, 'wb') as out:
        out.write(down.read())
    try:
        unpack_archive(filename)
    except ReadError:
        pass
    return True


def deploy(software):
    """
    According to system, install software.
    Return False if failed
    """
    log.info('Try to install {}.'.format(software))
    log.warning('Please consider to install it following official '
                'instruction to get a CLEAN system.')
    sys = system()
    # url dict
    blast_url = ('ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.8.1/'
                 'ncbi-blast-2.8.1+')
    iqtree_url = ('https://github.com/Cibiv/IQ-TREE/releases/download/v1.6.9/'
                  'iqtree-1.6.9')
    mafft_url = 'https://mafft.cbrc.jp/alignment/software/mafft'
    # windows blast path not sure
    urls = {'Linux':
                {'BLAST': {'url': blast_url+'-x64-linux.tar.gz',
                           'path': abspath('ncbi-blast-2.8.1+'+sep+'bin')},
                 'IQTREE': {'url': iqtree_url+'-Linux.tar.gz',
                            'path': abspath('iqtree-1.6.9-Linux'+sep+'bin')},
                 'MAFFT': {'url': mafft_url+'-7.407-linux.tgz',
                           'path': abspath('mafft-linux64')}},
            'Darwin':
                {'BLAST': {'url': blast_url+'-x64-macosx.tar.gz',
                           'path': abspath('ncbi-blast-2.8.1+'+sep+'bin')},
                 'IQTREE': {'url': iqtree_url+'-MacOSX.zip',
                            'path': abspath('iqtree-1.6.9-MacOSX'+sep+'bin')},
                 'MAFFT': {'url': mafft_url+'-7.407-mac.zip',
                           'path': abspath('mafft-mac')}},
            'Windows':
                {'BLAST': {'url': blast_url+'-win64.exe',
                           'path': abspath('cli')},
                 'IQTREE': {'url': iqtree_url+'-Windows.zip',
                            'path': abspath('iqtree-1.6.9-Windows'+sep+'bin')},
                 'MAFFT': {'url': mafft_url+'-7.409-win64-signed.zip',
                           'path': abspath('mafft-win')}}}
    url = urls[sys][software]['url']
    # down
    if sys == 'Windows':
        if not download_software(url):
            return None
        if software == 'BLAST':
            log.info('BLAST on Windows has to been installed manually.')
            log.info('Please configure PATH variable correctly.')
            run('ncbi-blast-2.8.1+-win64.exe', shell=True)
    elif sys == 'Linux':
        ok = False
        for pack_mgr in ('apt', 'dnf', 'yum', 'pkg'):
            r = run('sudo {} install ncbi-blast+ iqtree mafft'.format(
                pack_mgr), shell=True)
            if r.returncode == 0:
                ok = True
                break
        if not ok:
            log.info('Cannot install {} with package manager. Try to '
                     'download.'.format(software))
            if not download_software(url):
                return None
    elif sys == 'Darwin':
        with open(devnull, 'w', encoding='utf-8') as f:
            r = run('brew --help', shell=True, stdout=f, stderr=f)
        if r.returncode == 0:
            run('brew install blast mafft brewsci/science/iqtree', shell=True)
        else:
            log.warning('Cannot find Homebrew.')
            if not download_software(url):
                return None
            # after unzip, file lost executable flag on mac system
            run('chmod +x {}'.format(join_path(urls[sys][software]['path'],
                                               '*')), shell=True)
            if software == 'MAFFT':
                run('chmod +x {}'.format(join_path(urls[sys]['MAFFT']['path'],
                                                   'mafftdir', 'bin', '*')),
                    shell=True)
                run('chmod +x {}'.format(join_path(urls[sys]['MAFFT']['path'],
                                                   'mafftdir', 'libexec',
                                                   '*')), shell=True)
    # windows can omit .bat, linux cannot
    if software == 'MAFFT' and sys != 'Windows':
        rename(join_path(urls[sys]['MAFFT']['path'], 'mafft.bat'),
               join_path(urls[sys]['MAFFT']['path'], 'mafft'))
    return abspath(urls[sys][software]['path'])