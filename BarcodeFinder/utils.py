#!/usr/bin/python3

import re
import logging
import platform

from functools import lru_cache
from pathlib import Path
from platform import system
from os import devnull as DEVNULL
from os.path import abspath, exists, pathsep, sep
from os.path import join as join_path
from os import (cpu_count, devnull, environ, pathsep, rename)
from urllib.error import HTTPError
from urllib.request import urlopen
from shutil import unpack_archive, ReadError
from subprocess import run

from Bio.Seq import Seq

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


class BlastResult:
    # slightly faster than namedtuple
    __slots = ('query_id', 'hit_id', 'query_seq', 'ident_num', 'mismatch_num',
               'bitscore_raw', 'query_start', 'query_end', 'hit_start',
               'hit_end')

    def __init__(self, line):
        record = line.strip().split('\t')
        self.query_id, self.hit_id, self.query_seq = record[0:3]
        (self.ident_num, self.mismatch_num, self.bitscore_raw,
         self.query_start, self.query_end, self.hit_start,
         self.hit_end) = [int(i) for i in record[3:]]

    def __repr__(self):
        fmt = ('query_id: {}\thit_id: {}\tbitscore_raw: {}\tquery_start:{}\t'
               'query_end: {}\thit_start: {}\thit_end: {}')
        return fmt.format(self.query_id, self.hit_id, self.bitscore_raw,
                          self.query_start, self.query_end, self.hit_start,
                          self.hit_end)


def init_out(arg):
    """
    Initilize output folder.
    Args:
        arg(NameSpace): arguments
    Returns:
        arg(NameSpace): arguments
    """
    if arg.out is None:
        log.warning('Output folder was not set.')
        log.info('\tUse "Result" instead.')
        arg.out = Path().cwd().absolute() / 'Result'
    else:
        arg.out = Path(arg.out).absolute()
    if arg.out.exists():
        log.error(f'Output folder {arg.out} exists.')
        return None
    try:
        arg.out.mkdir()
        arg._gb = arg.out / 'GenBank'
        arg._gb.mkdir()
        arg._fasta = arg.out / 'Fasta'
        arg._fasta.mkdir()
        arg._divide = arg.out / 'Divide'
        arg._divide.mkdir()
        arg._expand = arg.out / 'Expanded_fasta'
        arg._expand.mkdir()
        arg._uniq = arg.out / 'Uniq'
        arg._uniq.mkdir()
        arg._align = arg.out / 'Alignment'
        arg._align.mkdir()
        arg._evaluate = arg.out /'Evaluate'
        arg._evaluate.mkdir()
        arg._tmp = arg.out / 'Temp'
        arg._tmp.mkdir()
    except Exception:
        log.warning('Folder exists.')
    return arg


@lru_cache(maxsize=None)
def gene_rename(old_name: str, genbank_format=False) -> (str, str):
    """
    Old doc:
        Different name of same gene will cause data to be splited to numerous
        files instead of one and some data may be dropped.

        For chloroplast genes, the author summarized various kinds of
        annotation error of gene name or synonyms and try to use regular
        expression to fix it.

        Ideally, use BLAST to re-annotate sequence is the best(and slow) way to
        find the correct name. This function only offers a "hotfix" which is
        enough.
    Rename plastid genes.
    May be dangerous.
    Will cache results.
    Args:
        old_name: old gene name
        genbank_format: use style like "trnH-GUG" or "trnHgug"
    Returns:
        new_name(str): new name, if fail, return old name
        gene_type(str): gene types, guessed from name
    """
    lower = old_name.lower()
    # (trna|trn(?=[b-z]))
    s = re.compile(r'(\d+\.?\d?)(s|rrn|rdna)')
    if lower.startswith('trn'):
        pattern = re.compile(r'([atcgu]{3})')
        prefix = 'trn'
        aa_letter = 'X'
        try:
            anticodon = Seq(re.search(pattern, lower[3:]).group(1))
        except Exception:
            return old_name, 'bad_name'
        # rna editing? trnI-CAU
        if anticodon == 'cau' and lower.startswith('trni'):
            aa_letter = 'I'
        # for trnfM-CAU
        elif lower.startswith('trnfm'):
            prefix = 'trnf'
            aa_letter = 'M'
        else:
            aa_letter = anticodon.reverse_complement().translate().upper()
            #anticodon = anticodon.transcribe()
        if genbank_format:
            new_name = f'{prefix}{aa_letter}-{anticodon.upper()}'
        else:
            new_name = f'{prefix}{aa_letter}{anticodon.lower()}'
        gene_type = 'tRNA'
    elif lower.startswith('rrn'):
        pattern = re.compile(r'(\d+\.?\d?)')
        try:
            number = re.search(pattern, lower).group(1)
        except Exception:
            return old_name, 'bad_name'
        new_name = 'rrn{}'.format(number)
        gene_type = 'rRNA'
    elif re.search(s, lower) is not None:
        new_name = 'rrn{}'.format(re.search(s, lower).group(1))
        gene_type = 'rRNA'
    else:
        pattern = re.compile(r'[^a-z]*'
                             '(?P<gene>[a-z]+)'
                             '[^a-z0-9]*'
                             '(?P<suffix>[a-z]|[0-9]+)')
        match = re.search(pattern, lower)
        try:
            gene = match.group('gene')
            suffix = match.group('suffix')
        except Exception:
            return old_name, 'bad_name'
        new_name = '{}{}'.format(gene, suffix.upper())
        # capitalize last letter
        if len(new_name) > 3:
            s = list(new_name)
            if s[-1].isalpha():
                new_name = '{}{}'.format(
                    ''.join(s[:-1]), ''.join(s[-1]).upper())
        gene_type = 'normal'
    if len(lower) >= 15:
        gene_type = 'suspicious_name'
    return new_name, gene_type


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


def accessible(name: Path, type_: str) -> bool:
    """
    Check given path is accessible or not.
    Given path does not exist.
    Args:
        name(Path): folder or file, absolute path
        type_(str): 'folder' or 'file'
    Return:
        ok(bool): accessible or not
    """
    p = Path(name)
    if type_ == 'folder':
        try:
            p.mkdir()
            p.rmdir()
            ok = True
        except PermissionError:
            ok = False
    elif type_ == 'file':
        try:
            p.touch()
            p.unlink()
            ok = True
        except PermissionError:
            ok = False
    else:
        log.critical(f'Illegal type: {type_}')
        ok = False
    return ok


def test_cmd(program, option='-version'):
    """
    Test given program and option is ok to run or not.
    Args:
        program(Path or str): program path, could be relative path if it can
        be found in $PATH or %PATH%
        option(str): option for program, usually use "-v" to show version to
        test the program
    Return:
        success(bool): success or not
    """
    test = run(f'{program} {option}', shell=True, stdout=DEVNULL,
               stderr=DEVNULL)
    success = True if test.returncode == 0 else False
    return success


def get_third_party():
    """
    Get third_party folder.
    If do not exist, create it.
    If cannot access, report.
    Return:
        success(bool): ok or not
        third_party(Path): absolute path of third_party folder
    """
    third_party = Path().home().absolute() / '.barcodefinder'
    success = False
    if not third_party.exists():
        log.debug(f'Create folder {third_party}')
        try:
            third_party.mkdir()
        except Exception:
            log.critical(f'Failed to create {third_party}.'
                         'Please contact the administrator.')
            return success, third_party
    if not accessible(third_party/'test', 'file'):
        log.critical(f'Failed to access {third_party}.'
                     f'Please contact the administrator.')
        return success, third_party
    success = True
    return success, third_party


def get_blast(third_party=None):
    """
    Get BLAST location.
    If BLAST was found, assume makeblastdb is found, too.
    If not found, download it.
    Args:
        third_party(Path or None): path for install
    Return:
        ok(bool): success or not
        blast(str): blast path
    """
    if third_party is None:
        third_party_ok, third_party = get_third_party()
        if not third_party_ok:
            return third_party_ok, ''
    home_blast = third_party / 'ncbi-blast-2.10.0+' / 'bin' / 'blastn'
    # in Windows, ".exe" can be omitted
    # win_home_blast = home_blast.with_name('blastn.exe')
    ok = False
    # older than 2.8.1 is buggy
    url = ('ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.0/'
           'ncbi-blast-2.10.0+')
    urls = {'Linux': url+'-x64-linux.tar.gz',
            'Darwin': url+'-x64-macosx.tar.gz',
            'Windows': url+'-x64-win64.tar.gz'}
    blast = 'blastn'
    if test_cmd(blast)
        ok = True
        return ok, blast
    if test_cmd(home_blast)
        ok = True
        return ok, str(home_blast)
    log.warning('Cannot find NCBI BLAST, try to install.')
    log.info('According to Internet speed, may be slow.')
    try:
        # 50kb/10s=5kb/s, enough for test
        _ = urlopen('https://www.ncbi.nlm.nih.gov', timeout=10)
    except Exception:
        log.critical('Cannot connect to NCBI.')
        log.critical('Please check your Internet connection.')
        return ok, ''
    try:
        # file is 86-222mb, 222mb/3600s=60kb/s, consider it's ok for users all
        # over the world
        down = urlopen(urls[platform.system()], timeout=3600)
    except Exception:
        log.critical('Cannot download BLAST.')
        log.critical('Please manually download it from'
                     'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/')
        return ok, ''
    down_file = third_party / 'BLAST_2.10.0.tar.gz'
    with open(down_file, 'wb') as out:
        out.write(down.read())
    try:
        unpack_archive(down_file, third_party)
    except Exception:
        log.critical('The file is damaged.')
        log.critical('Please check your connection.')
        return ok, ''
    assert test_cmd(home_blast, '-version')
    ok = True
    return ok, str(home_blast)


def get_all_third_party():
    """
    Use three threads to speed up.
    """
    log.info('Try to locate or install all third-party software.')
    third_party_ok, third_party = get_third_party()
    if not third_party_ok:
        return -1
    iqtree = Thread(target=get_iqtree, args=(third_party,), daemon=True)
    novoplasty = Thread(target=get_novoplasty, args=(third_party,),
                        daemon=True)
    blast = Thread(target=get_blast, args=(third_party,), daemon=True)
    perl.start()
    novoplasty.start()
    blast.start()
    perl.join()
    novoplasty.join()
    blast.join()
    log.info('Finished.')
    return


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


def parse_blast_tab(filename):
    """
    Parse BLAST result (tab format).
    """
    query = []
    with open(filename, 'r', encoding='utf-8') as raw:
        for line in raw:
            if line.startswith('# BLAST'):
                yield query
                query = []
            elif line.startswith('#'):
                pass
            else:
                query.append(BlastResult(line))



