#!/usr/bin/python3

import re
import logging
import platform
import subprocess

from collections import Iterable
from functools import lru_cache
from queue import Queue
from threading import Thread
from pathlib import Path
from urllib.request import urlopen
from shutil import unpack_archive

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


def arg_to_str(arg) ->str:
    s = ''
    arg_dict = vars(arg)
    for key, value in arg_dict.items():
        if isinstance(value, str):
            pass
        elif isinstance(value, Iterable):
            value = ' '.join(value)
        elif value is None:
            continue
        elif isinstance(value, bool):
            if value:
                value = ''
            else:
                # Assume all bool option using action='store_true'
                continue
        s += f' -{key} {value}'
    return s


def add_file_log(arg):
    """
    Add file handler if not exist.
    Note that "log" is a global variable.
    """
    has_file_hdl = any([type(i)==logging.FileHandler for i in log.handlers])
    if has_file_hdl:
        pass
    else:
        log_file = arg.out / 'Log.txt'
        log_file_handler = logging.FileHandler(log_file, mode='a')
        log_file_handler.setLevel(logging.INFO)
        log_file_handler.setFormatter(
            logging.Formatter( fmt=FMT, datefmt=DATEFMT))
        log.addHandler(log_file_handler)
        log.debug('Add file handler.')
    return


def move(source: Path, dest, copy=False):
    """
    Move source to dest and return dest.
    If set "copy", copy source to dest instead of move.
    Because Path.rename could not move file across different filesystem or
    drive, have to use copy and delete to implement "move".
    Warning:
        This function does not check whether dest exists or not.
    Args:
        source(Path): old path
        dest(Path or str): new path
        copy(bool): copy or move
    Return:
        dest(Path): new path
    """
    source = Path(source).absolute()
    dest = Path(dest).absolute()
    # avoid useless copy
    # Path.samefile may raise FileNotFoundError
    if source == dest:
        log.debug(f'{source} and {dest} are same.')
    else:
        # read_bytes/write_bytes includes open, read/write and close steps
        dest.write_bytes(source.read_bytes())
        if not copy:
            source.unlink()
    return dest


def init_out(arg):
    """
    Initilize output folder.
    Args:
        arg(NameSpace): arguments
    Returns:
        arg(NameSpace): arguments
    """
    if not hasattr(arg, 'out') or arg.out is None:
        log.warning('Output folder was not set.')
        log.info('\tUse "Result" instead.')
        arg.out = Path().cwd().absolute() / 'Result'
    else:
        arg.out = Path(arg.out).absolute()
    if arg.out.exists():
        from BarcodeFinder import global_vars
        if not global_vars.global_dict.get('out_inited', False):
            log.error(f'Output folder {arg.out} exists.')
            arg.out = arg.out.parent / (arg.out.name+'_')
            log.info(f'Use {arg.out} instead.')
            if arg.out.exists():
                log.error(f'{arg.out} exists, too!')
                raise SystemExit(-1)
        else:
            pass
    arg._gb = arg.out / 'GenBank'
    arg._fasta = arg.out / 'Fasta'
    arg._divide = arg.out / 'Divide'
    arg._expand = arg.out / 'Expanded_fasta'
    arg._unique = arg.out / 'Unique'
    arg._align = arg.out / 'Alignment'
    arg._evaluate = arg.out /'Evaluate'
    arg._primer = arg.out / 'Primer'
    arg._tmp = arg.out / 'Temp'
    try:
        arg.out.mkdir()
        arg._gb.mkdir()
        arg._fasta.mkdir()
        arg._divide.mkdir()
        arg._expand.mkdir()
        arg._unique.mkdir()
        arg._align.mkdir()
        arg._evaluate.mkdir()
        arg._primer.mkdir()
        arg._tmp.mkdir()
    except Exception:
        log.debug('Folder exists.')
    out_inited = True
    return arg


def clean_tmp(filename: Path):
    if filename.is_dir():
        for i in filename.glob('*'):
            i.unlink()
    elif filename.is_file():
        for i in filename.parent.glob(f'{filename.name}*'):
            i.unlink()
    return


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


def test_cmd(program, option='-version') -> bool:
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
    test = subprocess.run(f'{program} {option}', shell=True,
                          stdout=subprocess.DEVNULL,
                          stderr=subprocess.DEVNULL)
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


def get_blast(third_party=None, result=None) -> (bool, str):
    """
    Get BLAST location.
    If BLAST was found, assume makeblastdb is found, too.
    If not found, download it.
    Args:
        third_party(Path or None): path for install
        result(Queue): return values
    Return:
        ok(bool): success or not
        blast(str): blast path
    """
    if third_party is None:
        third_party_ok, third_party = get_third_party()
        if not third_party_ok:
            return third_party_ok, ''
    blast = 'blastn'
    home_blast = third_party / 'ncbi-blast-2.11.0+' / 'bin' / blast
    # in Windows, ".exe" can be omitted
    # win_home_blast = home_blast.with_name('blastn.exe')
    ok = False
    # older than 2.8.1 is buggy
    url = ('ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/'
           'ncbi-blast-2.11.0+')
    urls = {'Linux': url+'-x64-linux.tar.gz',
            'Darwin': url+'-x64-macosx.tar.gz',
            'Windows': url+'-x64-win64.tar.gz'}
    if test_cmd(blast):
        ok = True
        home_blast = str(blast)
    elif test_cmd(home_blast):
        ok = True
        home_blast = str(home_blast)
    else:
        while True:
            log.warning('Cannot find NCBI BLAST, try to install.')
            log.info('According to Internet speed, may be slow.')
            try:
                # 50kb/10s=5kb/s, enough for test
                _ = urlopen('https://www.ncbi.nlm.nih.gov', timeout=10)
            except Exception:
                log.critical('Cannot connect to NCBI.')
                log.critical('Please check your Internet connection.')
                break
            try:
                # file is 86-222mb
                down = urlopen(urls[platform.system()], timeout=10)
            except Exception:
                log.critical('Cannot download BLAST.')
                log.critical('Please manually download it from'
                             'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/')
                break
            down_file = third_party / 'BLAST_2.11.0.tar.gz'
            with open(down_file, 'wb') as out:
                out.write(down.read())
            try:
                unpack_archive(down_file, third_party)
            except Exception:
                log.critical('The file is damaged.')
                log.critical('Please check your connection.')
                break
            assert test_cmd(home_blast, '-version')
            ok = True
            break
    if result is not None:
        result.put(('BLAST', ok))
    return ok, str(home_blast)


def get_iqtree(third_party=None, result=None) -> (bool, str):
    """
    Get iqtree location.
    If not found, download it.
    Args:
        third_party(Path or None): path for install
        result(Queue): return values
    Return:
        ok(bool): success or not
        iqtree(str): blast path
    """
    if third_party is None:
        third_party_ok, third_party = get_third_party()
        if not third_party_ok:
            return third_party_ok, ''
    iqtree = 'iqtree2'
    # in Windows, ".exe" can be omitted
    ok = False
    url = 'https://github.com/Cibiv/IQ-TREE/releases/download/v2.0.6/'
    # platform: (filename, folder)
    fileinfo = {'Linux': ('iqtree-2.0.6-Linux.tar.gz', 'iqtree-2.0.6-Linux'),
                'Darwin': ('iqtree-2.0.6-MacOSX.zip', 'iqtree-2.0.6-MacOSX'),
                'Windows': ('iqtree-2.0.6-Windows.zip',
                            'iqtree-2.0.6-Windows')}
    system = platform.system()
    filename = fileinfo[system][0]
    home_iqtree = third_party / fileinfo[system][1] / 'bin' / iqtree
    if test_cmd(iqtree):
        ok = True
        home_iqtree = str(iqtree)
    elif test_cmd(home_iqtree):
        ok = True
    else:
        while True:
            log.warning('Cannot find iqtree, try to install.')
            log.info('According to Internet speed, may be slow.')
            try:
                # 50kb/10s=5kb/s, enough for test
                _ = urlopen('https://github.com', timeout=10)
            except Exception:
                log.critical('Cannot connect to github.com')
                log.critical('Please check your Internet connection.')
                break
            try:
                # file is ~10mb
                down = urlopen(f'{url}{filename}', timeout=10)
            except Exception:
                log.critical('Cannot download iqtree.')
                log.critical('Please manually download it from '
                             'https://github.com/Cibiv/IQ-TREE/')
                break
            down_file = third_party / fileinfo[system][0]
            with open(down_file, 'wb') as out:
                out.write(down.read())
            try:
                #unpack_archive(down_file, third_party/fileinfo[system][1])
                unpack_archive(down_file, third_party)
            except Exception:
                log.critical('The file is damaged.')
                log.critical('Please check your connection.')
                break
            assert test_cmd(home_iqtree, '-version')
            ok = True
            break
    if result is not None:
        result.put(('IQTREE', ok))
    return ok, str(home_iqtree)


def get_mafft(third_party=None, result=None) -> (bool, str):
    """
    Get iqtree location.
    If not found, download it.
    Args:
        third_party(Path or None): path for install
        result(Queue): return values
    Return:
        ok(bool): success or not
        iqtree(str): blast path
    """
    if third_party is None:
        third_party_ok, third_party = get_third_party()
        if not third_party_ok:
            return third_party_ok, ''
    system = platform.system()
    if system == 'Windows':
        mafft = 'mafft.bat'
    else:
        mafft = 'mafft'
    # in Windows, ".exe" can be omitted
    # win_home_blast = home_blast.with_name('blastn.exe')
    ok = False
    url = 'https://mafft.cbrc.jp/alignment/software/'
    fileinfo = {'Linux': ('mafft-7.475-linux.tgz', 'mafft-linux64'),
                'Darwin': ('mafft-7.475-mac.zip', 'mafft-mac'),
                'Windows': ('mafft-7.475-win64-signed.zip', 'mafft-win')}
    home_mafft = third_party / fileinfo[system][1] / mafft
    if test_cmd(mafft, '--version'):
        ok = True
        home_mafft = str(mafft)
    elif test_cmd(home_mafft, '--version'):
        ok = True
    else:
        while True:
            log.warning('Cannot find mafft, try to install.')
            log.info('According to Internet speed, may be slow.')
            try:
                # 50kb/10s=5kb/s, enough for test
                _ = urlopen('https://mafft.cbrc.jp', timeout=10)
            except Exception:
                log.critical('Cannot connect to mafft.cbrc.jp')
                log.critical('Please check your Internet connection.')
                break
            try:
                # file is ~10mb
                if system != 'Darwin':
                    down = urlopen(url+fileinfo[system][0], timeout=10)
                else:
                    down = urlopen(url+fileinfo[system][0]+'?signed', timeout=10)
            except Exception:
                log.critical('Cannot download mafft.')
                log.critical(f'Please manually download it from {url}')
                break
            down_file = third_party / fileinfo[system][0]
            with open(down_file, 'wb') as out:
                out.write(down.read())
            try:
                #unpack_archive(down_file, third_party/fileinfo[system][1])
                unpack_archive(down_file, third_party)
            except Exception:
                log.critical('The file is damaged.')
                log.critical('Please check your connection.')
                break
            assert test_cmd(home_mafft, '--version')
            ok = True
            break
    if result is not None:
        result.put(('MAFFT', ok))
    return ok, str(home_mafft)


def get_all_third_party() -> bool:
    """
    Use three threads to speed up.
    """
    log.info('Try to locate or install all third-party software.')
    third_party_ok, third_party = get_third_party()
    if not third_party_ok:
        return False
    result_queue = Queue()
    iqtree = Thread(target=get_iqtree, args=(third_party, result_queue),
                    daemon=True)
    mafft = Thread(target=get_mafft, args=(third_party, result_queue),
                   daemon=True)
    blast = Thread(target=get_blast, args=(third_party, result_queue),
                   daemon=True)
    iqtree.start()
    mafft.start()
    blast.start()
    iqtree.join()
    mafft.join()
    blast.join()
    while not result_queue.empty():
        name, ok = result_queue.get()
        if ok:
            log.info(f'Got {name}.')
        else:
            log.error(f'Failed to got {name}.')
            return False
    return True


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
    pass


def codon_usage(alignment):
    pass


def gap_analyze(gap_alignment):
    """

    Args:
        gap_alignment: np.array

    Returns:

    """
    pass
