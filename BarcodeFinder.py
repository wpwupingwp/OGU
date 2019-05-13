#!/usr/bin/python3

import argparse
import logging
import json
import re
from collections import defaultdict
from glob import glob
from itertools import product as cartesian_product
from io import StringIO
from os import (cpu_count, devnull, environ, mkdir, pathsep, remove, rename,
                sep)
from os.path import abspath, basename, exists, splitext
from os.path import join as join_path
from platform import system
from random import choice
from shutil import unpack_archive, ReadError
from subprocess import run
from urllib.error import HTTPError
from urllib.request import urlopen

import coloredlogs
import numpy as np
from primer3 import calcTm, calcHairpinTm, calcHomodimerTm, calcHeterodimerTm
from Bio import Entrez, Phylo, SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline as Blast
from Bio.Data.IUPACData import ambiguous_dna_values as ambiguous_data
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

from matplotlib import use as mpl_use
if environ.get('DISPLAY', '') == '':
    mpl_use('Agg')
from matplotlib import pyplot as plt
from matplotlib import rcParams

rcParams['axes.labelsize'] = 16
rcParams['axes.linewidth'] = 1.5
rcParams['axes.titlesize'] = 25
rcParams['font.size'] = 16
rcParams['lines.linewidth'] = 1.5


# define logger
FMT = '%(asctime)s %(levelname)-8s %(message)s'
DATEFMT = '%I:%M:%S'
TEMP_LOG = 'Temp.log'
LOG_FMT = logging.Formatter(fmt=FMT, datefmt=DATEFMT)
log = logging.getLogger(__name__)
coloredlogs.install(level=logging.INFO, fmt=FMT, datefmt=DATEFMT)
# temporary filehandler before arg.out is set
log_tmp = logging.FileHandler(TEMP_LOG, 'w')
log_tmp.setLevel(logging.INFO)
log_tmp.setFormatter(LOG_FMT)
log.addHandler(log_tmp)


class PrimerWithInfo(SeqRecord):
    # inherit from Bio.SeqRecord.SeqRecord
    def __init__(self, seq='', quality=None, start=0, coverage=0,
                 avg_bitscore=0, mid_loc=None, avg_mismatch=0, detail=0,
                 is_reverse_complement=False):
        # store str
        super().__init__(Seq(seq.upper()))
        self.sequence = str(self.seq)

        # primer3.setGlobals seems have no effect on calcTm, use
        # calc_ambiguous_seq
        self.quality = self.letter_annotations['solexa_quality'] = quality
        self.start = self.annotations['start'] = start
        self.end = self.annotations['end'] = start + self.__len__() - 1
        self.coverage = self.annotations['coverage'] = coverage
        self.avg_bitscore = self.annotations['avg_bitscore'] = avg_bitscore
        self.mid_loc = self.annotations['mid_loc'] = mid_loc
        self.avg_mismatch = self.annotations['avg_mismatch'] = avg_mismatch
        self.detail = self.annotations['detail'] = detail
        self.is_reverse_complement = self.annotations['is_reverse_complement'
                                                      ] = False
        self.description = self.annotations['description'] = ''
        self.avg_mid_loc = 0
        self.hairpin_tm = 0
        self.homodimer_tm = 0
        self.tm = 0
        self.update_id()

    def __getitem__(self, i):
        # part of attribution do not change, others were reset
        if isinstance(i, int):
            i = slice(i, i + 1)
        if isinstance(i, slice):
            answer = PrimerWithInfo(seq=str(self.seq[i]),
                                    quality=self.quality[i])
            answer.annotations = dict(self.annotations.items())
            return answer
        else:
            log.exception('Bad index.')
            raise

    def reverse_complement(self):
        table = str.maketrans('ACGTMRWSYKVHDBXN', 'TGCAKYWSRMBDHVXN')
        new_seq = str.translate(self.sequence, table)[::-1]
        new_quality = self.quality[::-1]
        # try to simplify??
        return PrimerWithInfo(seq=new_seq, quality=new_quality,
                              start=self.start, coverage=self.coverage,
                              avg_bitscore=self.avg_bitscore,
                              mid_loc=self.mid_loc,
                              is_reverse_complement=True,
                              detail=self.detail)

    def update_id(self):
        self.end = self.annotations['end'] = self.start + self.__len__() - 1
        if self.mid_loc is not None and len(self.mid_loc) != 0:
            self.avg_mid_loc = int(average(list(self.mid_loc.values())))
        self.id = ('AvgMidLocation({:.0f})-Tm({:.2f})-Coverage({:.2%})-'
                   'AvgBitScore({:.2f})-Start({})-End({})'.format(
                    self.avg_mid_loc, self.tm, self.coverage,
                    self.avg_bitscore, self.start, self.end))


class Pair:
    # save memory
    __slots__ = ['left', 'right', 'delta_tm', 'coverage', 'start', 'end',
                 'resolution', 'tree_value', 'avg_terminal_len', 'entropy',
                 'have_heterodimer', 'heterodimer_tm', 'pi', 'score',
                 'length', 'gap_ratio']

    def __init__(self, left, right, alignment):
        rows, columns = alignment.shape
        self.left = left
        self.right = right
        self.delta_tm = abs(self.left.tm - self.right.tm)
        # get accurate length
        a = len(self.left) / 2
        b = len(self.right) / 2
        common = left.mid_loc.keys() & right.mid_loc.keys()
        lengths = [[key, ((right.mid_loc[key] - b) - (left.mid_loc[key] + a))
                    ] for key in common]
        lengths = {i[0]: int(i[1]) for i in lengths if i[1] > 0}
        self.length = lengths
        self.right.coverage = len(self.right.mid_loc) / rows
        self.coverage = len(common) / rows
        # pairs use mid_loc from BLAST as start/end
        self.start = self.left.avg_mid_loc
        self.end = self.right.avg_mid_loc
        self.have_heterodimer = False
        self.heterodimer_tm = 0.0
        self.resolution = 0.0
        self.tree_value = 0.0
        self.avg_terminal_len = 0.0
        self.entropy = 0.0
        self.pi = 0.0
        self.gap_ratio = 0.0
        self.score = self.get_score()

    def __repr__(self):
        return (
            'Pair(score={:.2f}, product={:.0f}, start={}, end={}, left={}, '
            'right={}, observerd_resolution={:.2%}, coverage={:.2%},'
            'delta_tm={:.2f}, have_heterodimer={})'.format(
                self.score, average(list(self.length.values())), self.start,
                self.end, self.left.seq, self.right.seq, self.resolution,
                self.coverage, self.delta_tm, self.have_heterodimer))

    def get_score(self):
        # calculate score of given primer pairs. Suggestion only
        # use score to filter primer pairs
        return (average(list(self.length.values())) * 0.5
                + self.coverage * 200
                + len(self.left) * 10
                + len(self.right) * 10
                + self.resolution * 100
                + self.tree_value * 100 + self.entropy * 5
                - int(self.have_heterodimer) * 10
                - self.delta_tm * 5 - self.left.avg_mismatch * 10
                - self.right.avg_mismatch * 10)

    def add_info(self, alignment):
        # put attributes that need heavy computation here for the final primer
        # pairs in order to save CPU time
        if not self.right.is_reverse_complement:
            self.right = self.right.reverse_complement()
        # include end base, use alignment loc for slice
        (self.gap_ratio, self.resolution, self.entropy, self.pi,
         self.tree_value, self.avg_terminal_len) = get_resolution(
             alignment, self.left.start, self.right.end + 1)
        self.heterodimer_tm = calc_ambiguous_seq(calcHeterodimerTm,
                                                 self.left.seq, self.right.seq)
        if max(self.heterodimer_tm, self.left.tm,
               self.right.tm) == self.heterodimer_tm:
            self.have_heterodimer = True
        else:
            self.have_heterodimer = False
        self.get_score()
        self.left.update_id()
        self.right.update_id()
        return self


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


def parse_args():
    """
    Parse args and store some global/temporary values.
    """
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=main.__doc__)
    general = arg.add_argument_group('General')
    general.add_argument('-aln', help='aligned fasta files to analyze')
    general.add_argument('-fasta', help='unaligned fasta format data to add')
    general.add_argument('-gb', help='genbank files')
    general.add_argument('-stop', type=int, choices=(1, 2),
                         help=('Stop after which step:'
                               '\t1. Download and pre-process;'
                               '\t2. Analyze variance;'))
    general.add_argument('-out', help='output directory')
    genbank = arg.add_argument_group('Genbank')
    genbank.add_argument('-email', help='email address for querying Genbank')
    genbank.add_argument('-exclude', help='exclude option')
    genbank.add_argument('-gene', type=str, help='gene name')
    genbank.add_argument('-group',
                         choices=('animals', 'plants', 'fungi', 'protists',
                                  'bacteria', 'archaea', 'viruses'),
                         help='Species kind')
    genbank.add_argument('-min_len', default=100, type=int,
                         help='minium length')
    genbank.add_argument('-max_len', default=10000, type=int,
                         help='maximum length')
    genbank.add_argument('-molecular', choices=('DNA', 'RNA'),
                         help='molecular type')
    genbank.add_argument('-og', '-organelle', dest='organelle',
                         choices=('mt', 'mitochondrion', 'cp',
                                  'chloroplast', 'pl', 'plastid'),
                         help='organelle type')
    genbank.add_argument('-query', nargs='*', help='query text')
    genbank.add_argument('-refseq', action='store_true',
                         help='Only search in RefSeq database')
    genbank.add_argument('-seq_n', default=None, type=int,
                         help='maximum number of records to download')
    genbank.add_argument('-taxon', help='Taxonomy name')
    pre = arg.add_argument_group('Preprocess')
    pre.add_argument('-expand', type=int, default=200,
                     help='expand length of upstream/downstream')
    pre.add_argument('-max_name_len', default=50,
                     help='maximum length of feature name')
    pre.add_argument('-max_seq_len', default=20000,
                     help='maximum length of feature sequence')
    pre.add_argument('-no_divide', action='store_true',
                     help='analyze whole sequence instead of divided fragment')
    pre.add_argument('-rename', action='store_true', help='try to rename gene')
    pre.add_argument('-uniq', choices=('longest', 'random', 'first', 'no'),
                     default='first',
                     help='method to remove redundant sequences')
    evaluate = arg.add_argument_group('Evaluate')
    evaluate.add_argument('-fast', action='store_true', default=False,
                          help='faster evaluate variance by omit tree_value'
                          'and terminal branch length')
    evaluate.add_argument('-step', type=int, default=50,
                          help='step length for sliding-window scan')
    primer = arg.add_argument_group('Primer')
    primer.add_argument('-a', dest='ambiguous_base_n', type=int, default=4,
                        help='number of ambiguous bases')
    primer.add_argument('-c', dest='coverage', type=float, default=0.6,
                        help='minium coverage of base and primer')
    primer.add_argument('-m', dest='mismatch', type=int, default=4,
                        help='maximum mismatch bases in primer')
    primer.add_argument('-pmin', dest='min_primer', type=int, default=18,
                        help='minimum primer length')
    primer.add_argument('-pmax', dest='max_primer', type=int, default=24,
                        help='maximum primer length')
    primer.add_argument('-r', dest='resolution', type=float, default=0.5,
                        help='minium resolution')
    primer.add_argument('-t', dest='top_n', type=int, default=1,
                        help='keep n primers for each high varient region')
    primer.add_argument('-tmin', dest='min_product', type=int, default=350,
                        help='minimum product length(include primer)')
    primer.add_argument('-tmax', dest='max_product', type=int, default=600,
                        help='maximum product length(include primer)')
    parsed = arg.parse_args()
    if parsed.fast:
        log.info('The "-fast" mode was opened. '
                 'Skip sliding-window scan with tree.')
    if parsed.group is not None:
        log.warning('The filters "group" was reported to return abnormal '
                    'records by Genbank. Please consider to use "-taxon" '
                    'instead.')
    # join nargs
    if parsed.query is not None:
        parsed.query = ' '.join(parsed.query)
    if parsed.refseq:
        log.info('Reset the limitation of sequence length for RefSeq.')
        parsed.min_len = None
        parsed.max_len = None
    if parsed.rename:
        log.warning('BarcodeFinder will try to rename genes by regular '
                    'expression.')
    if not any([parsed.query, parsed.taxon, parsed.group, parsed.gene,
                parsed.fasta, parsed.aln, parsed.gb, parsed.organelle]):
        return None
    if parsed.out is None:
        log.warning('Output folder was not set.')
        log.info('Use "Result" instead.')
        parsed.out = 'Result'
    parsed.by_gene_folder = join_path(parsed.out, 'by-gene')
    parsed.by_name_folder = join_path(parsed.out, 'by-name')
    # temporary filename, omit one parameters in many functions
    parsed.db_file = join_path(parsed.out, 'interleaved.fasta')
    parsed.out_file = ''
    # load option.json may cause chaos, remove
    return parsed


def average(x):
    """
    Safe average.
    """
    if len(x) == 0:
        return 0
    else:
        return sum(x) / len(x)


def safe(old):
    """
    Remove illegal character in file path or name.
    """
    return re.sub(r'\W', '_', old)


def clean_path(old, arg):
    """
    Join path if the file is not under by-gene or by-uniq to make working
    folder clean.
    """
    split = old.split(sep)
    if 'by-gene' not in split and 'by-name' not in split:
        return join_path(arg.by_name_folder, basename(old))
    else:
        return old


def calc_ambiguous_seq(func, seq, seq2=None):
    """
    Expand sequences with ambiguous bases to several clean sequences and apply
    func to every sequence.
    Return average value. Return 0 if len(seq) > 60 (from primer3)
    """
    # Seems primer3 only accept seqs shorter than 60 bp. Plus, too long seq
    # will cost too much memory.
    LEN_LIMIT = 60

    def _expand(seq):
        seq_list = []
        for base in seq:
            # replace illegal base with 'N'
            if base not in ambiguous_data:
                base = 'N'
            seq_list.append(ambiguous_data[base])
        seq_product = list(cartesian_product(*seq_list))
        seq_str = [''.join(i) for i in seq_product]
        return seq_str

    if len(seq) > LEN_LIMIT:
        log.warning('Ambiguous sequences is too long. Skip')
        return 0
    seq_str = _expand(seq)
    if seq2 is None:
        values = [func(i) for i in seq_str]
    else:
        if len(seq2) > LEN_LIMIT:
            log.warning('Ambiguous sequences is too long. Skip')
            return 0
        seq_str2 = _expand(seq2)
        products = cartesian_product(seq_str, seq_str2)
        values = [func(i[0], i[1]) for i in products]
    # primer3 will return negative values sometime
    values_positive = [max(0, i) for i in values]
    return average(values_positive)


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
        log.warn('Cannot download {}.'.format(filename))
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
                       'path': abspath('.')},
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


def get_query_string(arg):
    """
    Based on given options, generate query string from Genbank.
    """
    condition = []
    if arg.group is not None:
        condition.append('{}[filter]'.format(arg.group))
    if arg.query is not None:
        condition.append(arg.query)
    if arg.gene is not None:
        if ' ' in arg.gene:
            condition.append('"{}"[gene]'.format(arg.gene))
        else:
            condition.append('{}[gene]'.format(arg.gene))
    if arg.molecular is not None:
        d = {'DNA': 'biomol_genomic[PROP]',
             'RNA': 'biomol_mrna[PROP]'}
        condition.append(d[arg.molecular])
    if arg.taxon is not None:
        condition.append('{}[ORGANISM]'.format(arg.taxon))
    if arg.organelle is not None:
        if arg.organelle in ('mt', 'mitochondrion'):
            condition.append('{mitochondrion}[filter]')
        else:
            condition.append('(plastid[filter] OR chloroplast[filter])')
    if arg.refseq:
        condition.append('refseq[filter]')
    if (len(condition) > 0) and (arg.min_len is not None and arg.max_len is
                                 not None):
        condition.append('("{}"[SLEN] : "{}"[SLEN])'.format(
            arg.min_len, arg.max_len))
    if arg.exclude is not None:
        condition.append('NOT ({})'.format(arg.exclude))
    if not condition:
        return None
    else:
        string = ' AND '.join(condition)
        string = string.replace('AND NOT', 'NOT')
        return string


def download(arg, query):
    """
    Download records from Genbank.
    Because of connection to Genbank website is not stable (especially in
    Asia), it will retry if failed. Ctrl+C to break.
    """

    TOO_MUCH = 50000
    # 100 is enough?
    RETRY_MAX = 100
    log.info('Query:\t{}'.format(query))
    if arg.email is None:
        Entrez.email = 'guest@example.com'
        log.warning('You did not provide email address for using Entrez, '
                    'use {} instead.'.format(Entrez.email))
    else:
        Entrez.email = arg.email
    query_handle = Entrez.read(Entrez.esearch(db='nuccore', term=query,
                                              usehistory='y'))
    count = int(query_handle['Count'])

    if count == 0:
        log.warning('Got 0 record. Please check the query.')
        log.info('Abort download.')
        return None
    elif count > TOO_MUCH:
        log.warning('Got {} records. Please note that if you repeatedly '
                    'sending huge queries you may be blocked by NCBI.'.format(
                        count))
    else:
        log.info('Got {} records.'.format(count))
    if arg.seq_n is not None:
        if count > arg.seq_n:
            count = arg.seq_n
            log.info('Download {} records becase of "-seq_n".'.format(
                arg.seq_n))
    log.info('Downloading... Ctrl+C to quit.')
    name_words = []
    for i in (arg.group, arg.taxon, arg.organelle, arg.gene, arg.query):
        if i is not None:
            name_words.append(i)
    if len(name_words) != 0:
        name = safe('-'.join(name_words))
    else:
        name = 'sequence'
    name = safe(name)
    file_name = join_path(arg.out, name + '.gb')
    output = open(file_name, 'w', encoding='utf-8')
    ret_start = 0
    if count >= 1000:
        ret_max = 1000
    elif count >= 100:
        ret_max = 100
    else:
        ret_max = 10
    retry = 0
    while ret_start < count:
        log.info('{:d}--{:d}'.format(ret_start, ret_start + ret_max))
        try:
            data = Entrez.efetch(db='nuccore',
                                 webenv=query_handle['WebEnv'],
                                 query_key=query_handle['QueryKey'],
                                 rettype='gb',
                                 retmode='text',
                                 retstart=ret_start,
                                 retmax=ret_max)
            output.write(data.read())
        # just retry if connection failed
        except IOError:
            if retry <= RETRY_MAX:
                log.warning('Failed on download. Retrying...')
                retry += 1
                continue
            else:
                log.critical('Too much failure ({} times).'.format(RETRY_MAX))
                log.info('Abort download.')
                return None
        ret_start += ret_max
    log.info('Download finished.')
    json_file = join_path(arg.out, 'Query.json')
    with open(json_file, 'w', encoding='utf-8') as _:
        json.dump(query_handle, _, indent=4, sort_keys=True)
    log.info('The query info was dumped into {}'.format(json_file))
    return file_name


def clean_gb(gbfile):
    """
    Records in Genbank may be problematic. Check it before parse and skip
    abnormal records.
    """
    log.info('Check Genbank file to remove abnormal records.')

    def parse_gb(handle):
        record = []
        for line in handle:
            record.append(line)
            if line.startswith('//'):
                yield record
                record = []
            else:
                pass

    wrong = 0
    old_gb = open(gbfile, 'r')
    for record in parse_gb(old_gb):
        # StringIO is faster than write tmp file to disk and read
        tmp_gb = StringIO()
        for _ in record:
            tmp_gb.write(_)
        tmp_gb.seek(0)
        try:
            gb_record = SeqIO.read(tmp_gb, 'gb')
            yield gb_record
        except ValueError as e:
            log.critical('Found problematic record {}: {}'.format(
                record[0][:25], e.args[0]))
            wrong += 1
    tmp_gb.close()
    old_gb.close()
    if wrong != 0:
        log.info('Remove {} abnormal records.'.format(wrong))


def gene_rename(old_name):
    """
    Different name of same gene will cause data to be splited to numerous
    files instead of one and some data may be dropped.
    For chloroplast genes, the auther summarized various kinds of annotation
    error of gene name or synonyms and try to use regular expression to fix
    it.
    Ideally, use BLAST to re-annotate sequence is the best(and slow) way to
    find the correct name. This function only offers a "hotfix".
    """
    lower = old_name.lower()
    # (trna|trn(?=[b-z]))
    s = re.compile(r'(\d+\.?\d?)(s|rrn|rdna)')
    if lower.startswith('trn'):
        pattern = re.compile(r'([atcgu]{3})')
        search = re.search(pattern, lower)
        if search is not None:
            codon = Seq(search.group(1))
        else:
            return old_name, 'bad_name'
        try:
            new_name = 'trn{}{}'.format(codon.reverse_complement().translate(),
                                        codon.transcribe())
        except ValueError:
            return old_name, 'bad_name'
        gene_type = 'tRNA'
    elif lower.startswith('rrn'):
        pattern = re.compile(r'(\d+\.?\d?)')
        search = re.search(pattern, lower)
        if search is not None:
            number = search.group(1)
        else:
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
        if match is not None:
            try:
                gene = match.group('gene')
                suffix = match.group('suffix')
            except ValueError:
                return old_name, 'bad_name'
        else:
            return old_name, 'bad_name'
        new_name = '{}{}'.format(gene, suffix.upper())
        # captitalize last letter
        if len(new_name) > 3:
            s = list(new_name)
            if s[-1].isalpha():
                new_name = '{}{}'.format(
                    ''.join(s[:-1]), ''.join(s[-1]).upper())
        gene_type = 'normal'
    # too long to be valid name
    if len(lower) >= 15:
        gene_type = 'suspicious_name'
    return new_name, gene_type


def write_seq(name, sequence_id, feature, whole_seq, path, arg):
    """
    Write fasta file.
    """

    def careful_extract(whole_seq):
        # illegal annotation may cause extraction failed
        try:
            sequence = feature.extract(whole_seq)
        except ValueError:
            sequence = ''
            log.warning('Cannot extract sequence of {} from {}.'.format(
                name, sequence_id))
        return sequence

    filename = join_path(path, name + '.fasta')
    sequence = careful_extract(whole_seq)
    with open(filename, 'a', encoding='utf-8') as handle:
        handle.write(sequence_id + '\n')
        handle.write(str(sequence) + '\n')
    if arg.expand != 0:
        if feature.location_operator == 'join':
            loc = feature.location.parts
            # ensure increasing order
            # parts do not have sort method
            loc.sort(key=lambda x: x.start)
            new_loc = sum([
                # avoid IndexError
                FeatureLocation(max(0, loc[0].start - arg.expand),
                                loc[0].end, loc[0].strand),
                *loc[1:-1],
                FeatureLocation(loc[-1].start,
                                min(len(whole_seq), loc[-1].end+arg.expand),
                                loc[-1].strand)])
            feature.location = new_loc
        feature.type = 'expand'
        sequence = careful_extract(whole_seq)
        filename2 = join_path(path, '{}.expand'.format(name))
        with open(filename2, 'a', encoding='utf-8') as handle:
            handle.write(sequence_id + '\n')
            handle.write(str(sequence) + '\n')
        return filename2
    return filename


def get_feature_name(feature, arg):
    """
    Get feature name and collect genes for extract spacer.
    Only handle gene, product, misc_feature, misc_RNA.
    """
    name = None
    if feature.type == 'gene':
        if 'gene' in feature.qualifiers:
            name = feature.qualifiers['gene'][0]
            if arg.rename:
                name = gene_rename(name)[0]
        elif 'product' in feature.qualifiers:
            name = feature.qualifiers['product'][0]
        elif 'locus_tag' in feature.qualifiers:
            name = feature.qualifiers['locus_tag'][0]
        else:
            log.warning('Cannot recognize annotation:\n{}'.format(feature))
    elif feature.type == 'misc_feature':
        if 'product' in feature.qualifiers:
            name = feature.qualifiers['product'][0]
        elif 'note' in feature.qualifiers:
            name = feature.qualifiers['note'][0]
        else:
            log.warning('Cannot recognize annotation:\n{}'.format(feature))
        if (name is not None) and ('intergenic_spacer' in name or 'IGS' in
                                   name):
            # 'IGS' in name) and len(name) < 100):
            name = name.replace('intergenic_spacer_region', 'IGS')
    elif feature.type == 'misc_RNA':
        if 'product' in feature.qualifiers:
            name = feature.qualifiers['product'][0]
        elif 'note' in feature.qualifiers:
            name = feature.qualifiers['note'][0]
        else:
            log.warning('Cannot recognize annotation:\n{}'.format(feature))
        # handle ITS
        if (name is not None) and 'internal_transcribed_spacer' in name:
            name = 'ITS'
    elif feature.type == 'rRNA':
        if 'product' in feature.qualifiers:
            name = feature.qualifiers['product'][0]
        elif 'note' in feature.qualifiers:
            name = feature.qualifiers['note'][0]
        else:
            log.warning('Cannot recognize annotation:\n{}'.format(feature))
    else:
        pass
    if name is not None:
        safe_name = safe(name)
    else:
        safe_name = None
    return safe_name, feature.type


def get_spacer(genes):
    """
    Given list of genes, extract spacers.
    """
    spacers = []
    # sorted according to sequence starting postion
    genes.sort(key=lambda x: int(x[1].location.start))
    for n, present in enumerate(genes[1:], 1):
        before = genes[n - 1]
        # use sort to handle complex location relationship of two fragments
        location = [before[1].location.start, before[1].location.end,
                    present[1].location.start, present[1].location.end]
        location.sort(key=lambda x: int(x))
        start, end = location[1:3]
        if before[1].location.strand == present[1].location.strand == -1:
            strand = -1
        else:
            strand = 1
        name = '_'.join([before[0], present[0]])
        spacer = SeqFeature(FeatureLocation(start, end), id=name,
                            type='spacer', strand=strand)
        spacers.append(spacer)
    return spacers


def divide(gbfile, arg):
    """
    Given genbank file, return divided fasta files.
    """
    # From NCBI Taxonomy Database, 2019.2.21 update
    # give all superkingdoms, kingdoms, phyla, classes because they do not
    # have uniform prefix
    log.info('Divide {} by annotation.'.format(gbfile))
    superkingdoms = {'Bacteria', 'Archaea', 'Eukaryota', 'Viruses', 'Viroids'}
    kingdoms = {'Fungi', 'Viridiplantae', 'Metazoa'}
    phyla_str = '''
    Bacteroidetes, Chlorobi, Cyanobacteria, Proteobacteria, Firmicutes,
    Deinococcus-Thermus, Xanthophyceae, Bacillariophyta, Phaeophyceae,
    Euglenida, Chlorophyta, Chytridiomycota, Ascomycota, Basidiomycota,
    Eustigmatophyceae, Apicomplexa, Microsporidia, Porifera, Cnidaria,
    Platyhelminthes, Nemertea, Nematoda, Annelida, Mollusca, Arthropoda,
    Brachiopoda, Echinodermata, Chordata, Rotifera, Ctenophora, Bryozoa,
    Rhombozoa, Hemichordata, Placozoa, Chaetognatha, Acanthocephala,
    Onychophora, Crenarchaeota, Euryarchaeota, Haplosporidia, Fusobacteria,
    Orthonectida, Nematomorpha, Gastrotricha, Priapulida, Streptophyta,
    Nitrospirae, Tardigrada, Entoprocta, Kinorhyncha, Candidatus Korarchaeota,
    Acidobacteria, Candidatus Marinimicrobia, Fibrobacteres, Gnathostomulida,
    Candidatus Bipolaricaulota, Candidatus Omnitrophica, Caldiserica,
    Candidatus Aminicenantes, Candidatus Atribacteria, Armatimonadetes,
    Candidatus Microgenomates, Dictyoglomi, Cycliophora, Candidatus
    Latescibacteria, Elusimicrobia, Verrucomicrobia, Bolidophyceae, Candidatus
    Saccharibacteria, Kiritimatiellaeota, Gemmatimonadetes, Candidatus
    Hydrogenedentes, Pinguiophyceae, Nanoarchaeota, Aquificae, Chloroflexi,
    Thermotogae, Deferribacteres, Chrysiogenetes, Thermodesulfobacteria,
    Actinobacteria, Planctomycetes, Spirochaetes, Chlamydiae, Candidatus
    Parcubacteria, Lentisphaerae, Candidatus Poribacteria, Loricifera,
    Candidatus Gracilibacteria, Picozoa, candidate division WWE3,
    Blastocladiomycota, Candidatus Cloacimonetes, Synergistetes, Tenericutes,
    Aurearenophyceae, candidate division NC10, Thaumarchaeota, Candidatus
    Aenigmarchaeota, Candidatus Diapherotrites, Chromerida, Candidatus
    Bathyarchaeota, Cryptomycota, Candidatus Calescamantes, Candidatus
    Aerophobetes, candidate division JL-ETNP-Z39, Candidatus Hydrothermae,
    Ignavibacteriae, candidate division WPS-1, candidate division WPS-2,
    Nitrospinae, Xenacoelomorpha, candidate division Zixibacteria, Candidatus
    Fervidibacteria, candidate division GAL15, Candidatus Geoarchaeota,
    candidate phylum NAG2, Candidatus Parvarchaeota, Colponemidia, Candidatus
    Berkelbacteria, candidate division CPR1, candidate division CPR2,
    candidate division CPR3, Candidatus Peregrinibacteria, Candidatus
    Lokiarchaeota, Candidatus Kapabacteria, candidate division WOR-3,
    candidate division KD3-62, Candidatus Thorarchaeota, Candidatus
    Fermentibacteria, Candidatus Rokubacteria, Candidatus Dadabacteria,
    Candidatus Curtissbacteria, Candidatus Daviesbacteria, Candidatus
    Levybacteria, Candidatus Gottesmanbacteria, Candidatus Shapirobacteria,
    Candidatus Woesebacteria, Candidatus Roizmanbacteria, Candidatus
    Pacebacteria, Candidatus Collierbacteria, Candidatus Beckwithbacteria,
    Candidatus Campbellbacteria, Candidatus Falkowbacteria, Candidatus
    Nomurabacteria, Candidatus Amesbacteria, Candidatus Magasanikbacteria,
    Candidatus Uhrbacteria, Candidatus Yanofskybacteria, Candidatus
    Kaiserbacteria, Candidatus Wolfebacteria, Candidatus Adlerbacteria,
    Candidatus Moranbacteria, Candidatus Giovannonibacteria, Candidatus
    Jorgensenbacteria, Candidatus Kuenenbacteria, Candidatus Azambacteria,
    Candidatus Melainabacteria, Candidatus Woesearchaeota, Candidatus
    Micrarchaeota, Candidatus Tectomicrobia, Candidatus Abawacabacteria,
    Candidatus Coatesbacteria, Candidatus Delongbacteria, Candidatus
    Doudnabacteria, Candidatus Edwardsbacteria, Candidatus Eisenbacteria,
    Candidatus Firestonebacteria, Candidatus Fischerbacteria, Candidatus
    Fraserbacteria, Candidatus Glassbacteria, Candidatus Handelsmanbacteria,
    Candidatus Lindowbacteria, Candidatus Margulisbacteria, Candidatus
    Raymondbacteria, Candidatus Riflebacteria, Candidatus Schekmanbacteria,
    Candidatus Wallbacteria, Candidatus Blackburnbacteria, Candidatus
    Woykebacteria, Candidatus Chisholmbacteria, Candidatus Andersenbacteria,
    Candidatus Brennerbacteria, Candidatus Buchananbacteria, Candidatus
    Colwellbacteria, Candidatus Harrisonbacteria, Candidatus Jacksonbacteria,
    Candidatus Kerfeldbacteria, Candidatus Komeilibacteria, Candidatus
    Liptonbacteria, Candidatus Lloydbacteria, Candidatus Nealsonbacteria,
    Candidatus Niyogibacteria, Candidatus Portnoybacteria, Candidatus
    Ryanbacteria, Candidatus Spechtbacteria, Candidatus Staskawiczbacteria,
    Candidatus Sungbacteria, Candidatus Tagabacteria, Candidatus
    Taylorbacteria, Candidatus Terrybacteria, Candidatus Veblenbacteria,
    Candidatus Vogelbacteria, Candidatus Wildermuthbacteria, Candidatus
    Yonathbacteria, Candidatus Zambryskibacteria, Candidatus Wirthbacteria,
    Candidatus Desantisbacteria, Rhodothermaeota, Candidatus Kryptonia,
    Mucoromycota, Zoopagomycota, Candidatus Verstraetearchaeota,
    Calditrichaeota, Candidatus Odinarchaeota, Candidatus Heimdallarchaeota,
    Balneolaeota, Candidatus Marsarchaeota, Candidatus Goldbacteria, candidate
    division AD3, candidate division FCPU426, Candidatus Abyssubacteria,
    Candidatus Aureabacteria, Abditibacteriota, Coprothermobacterota,
    Olpidiomycota'''
    phyla = {i.strip() for i in phyla_str.split(sep=',')}
    classes_str = '''
    Gammaproteobacteria, Actinobacteria, Bangiophyceae, Florideophyceae,
    Chrysophyceae, Dinophyceae, Cryptophyta, Chlorophyceae, Bryopsida,
    Cycadopsida, Gnetopsida, Liliopsida, Oomycetes, Saccharomycetes,
    Ustilaginomycetes, Heterolobosea, Colpodea, Litostomatea, Prostomatea,
    Nassophorea, Oligohymenophorea, Demospongiae, Hydrozoa, Anthozoa, Cubozoa,
    Scyphozoa, Trematoda, Cestoda, Pilidiophora, Enopla, Polychaeta,
    Gastropoda, Bivalvia, Cephalopoda, Polyplacophora, Branchiopoda,
    Ostracoda, Malacostraca, Merostomata, Arachnida, Chilopoda, Diplopoda,
    Asteroidea, Ophiuroidea, Echinoidea, Holothuroidea, Ascidiacea,
    Chondrichthyes, Amphibia, Aves, Monogononta, Gymnolaemata, Enteropneusta,
    Archiacanthocephala, Calcarea, Alphaproteobacteria, Betaproteobacteria,
    Deltaproteobacteria, Epsilonproteobacteria, Ginkgoopsida, Diplura,
    Protura, Collembola, Appendicularia, Thaliacea, Acoela, Mollicutes,
    Chloroflexia, Scaphopoda, Ulvophyceae, Pterobranchia, Karyorelictea,
    Spirotrichea, Coscinodiscophyceae, Bacillariophyceae, Fragilariophyceae,
    Synurophyceae, Phylactolaemata, Crinoidea, Labyrinthulomycetes,
    Pedinophyceae, Pelagophyceae, Phyllopharyngea, Monogenea,
    Glaucocystophyceae, Raphidophyceae, Dictyochophyceae, Catenulida,
    Mammalia, Clitellata, Eutardigrada, Hyphochytriomycetes, Gordioida,
    Bdelloidea, Palaeacanthocephala, Eoacanthocephala, Insecta,
    Nemertodermatida, Caudofoveata, Chrysomerophyceae, Pycnogonida,
    Coniferopsida, Hexactinellida, Symphyla, Entorrhizomycetes, Pauropoda,
    Acantharea, Polycystinea, Thermodesulfobacteria, Deferribacteres,
    Hexanauplia, Cephalocarida, Trebouxiophyceae, Homoscleromorpha,
    Phaeothamniophyceae, Remipedia, Acidimicrobiia, Rubrobacteria,
    Coriobacteriia, Seisonidea, Heterotardigrada, Bacilli,
    Mesostigmatophyceae, Stenolaemata, Andreaeopsida, Sphagnopsida,
    Polytrichopsida, Takakiopsida, Andreaeobryopsida, Lingulata, Craniata,
    Rhynchonellata, Flavobacteriia, Sphingobacteriia, Chrysiogenetes, Enoplea,
    Chromadorea, Ichthyosporea, Zygnemophyceae, Chlorokybophyceae,
    Klebsormidiophyceae, Spartobacteria, Tentaculata, Arthoniomycetes,
    Dothideomycetes, Eurotiomycetes, Lecanoromycetes, Leotiomycetes,
    Pezizomycetes, Sordariomycetes, Pneumocystidomycetes,
    Schizosaccharomycetes, Taphrinomycetes, Neolectomycetes, Tremellomycetes,
    Agaricomycetes, Placididea, Agaricostilbomycetes, Microbotryomycetes,
    Pucciniomycetes, Turbellaria, Methanonatronarchaeia, Polyacanthocephala,
    Thermoprotei, Methanobacteria, Methanococci, Halobacteria, Thermoplasmata,
    Thermococci, Archaeoglobi, Methanopyri, Actinopteri, Marchantiopsida,
    Jungermanniopsida, Clostridia, Aquificae, Thermotogae, Deinococci,
    Laboulbeniomycetes, Orbiliomycetes, Thermomicrobia, Nectonematoida,
    Chlorobia, Heterotrichea, Micrognathozoa, Bacteroidia, Dictyoglomia,
    Fusobacteriia, Verrucomicrobiae, Planctomycetia, Spirochaetia, Nitrospira,
    Chlamydiia, Fibrobacteria, Acidobacteriia, Glomeromycetes,
    Gemmatimonadetes, Methanomicrobia, Solenogastres, Actinophryidae,
    Polypodiopsida, Hypermastigia, Anaerolineae, Chitinivibrionia,
    Dehalococcoidia, Coleochaetophyceae, Charophyceae, Gloeobacteria,
    Lichinomycetes, Solibacteres, Katablepharidophyta, Monoplacophora,
    Ktedonobacteria, Leiosporocerotopsida, Anthocerotopsida, Haplomitriopsida,
    Oedipodiopsida, Tetraphidopsida, Opitutae, Aconoidasida, Wallemiomycetes,
    Atractiellomycetes, Classiculomycetes, Cryptomycocolacomycetes,
    Cystobasidiomycetes, Mixiomycetes, Stylonematophyceae, Candidatus
    Marinamargulisbacteria, Endomicrobia, Chytridiomycetes,
    Monoblepharidomycetes, Neocallimastigomycetes, Blastocladiomycetes,
    Exobasidiomycetes, Dacrymycetes, Thermolithobacteria, Caldilineae,
    Erysipelotrichia, Holophagae, Synchromophyceae, Zetaproteobacteria,
    Mediophyceae, Elusimicrobia, Synergistia, Armophorea, Rhodellophyceae,
    Geoglossomycetes, Phycisphaerae, Caldisericia, Cytophagia, Ignavibacteria,
    Nitriliruptoria, Negativicutes, Compsopogonophyceae, Mamiellophyceae,
    Armatimonadia, Nanohaloarchaea, Tritirachiomycetes, Archaeorhizomycetes,
    Chthonomonadetes, Xylonomycetes, Nephroselmidophyceae, Plagiopylea,
    Conoidasida, Nitrospinia, Coniocybomycetes, Oligosphaeria, Lentisphaeria,
    Picomonadea, Cladistia, Ardenticatenia, Basidiobolomycetes,
    Neozygitomycetes, Entomophthoromycetes, Sagittoidea, Thermoflexia,
    Thermoleophilia, Aphelidea, Lycopodiopsida, Chlorodendrophyceae,
    Moniliellomycetes, Malasseziomycetes, Oligoflexia, Blastocatellia,
    Nitrososphaeria, Globothalamea, Fimbriimonadia, Limnochordia, Candidatus
    Thalassoarchaea, Palaeonemertea, Geminibasidiomycetes, Tissierellia,
    Hadesarchaea, Spiculogloeomycetes, Candidatus Peribacteria, Longimicrobia,
    Acidithiobacillia, Acidobacteria subdivision 6, Candidatus
    Lambdaproteobacteria, Candidatus Muproteobacteria, Balneolia,
    Rhodothermia, Chitinophagia, Chitinispirillia, Candidatus Methanomethylia,
    Nuda, Kiritimatiellae, Staurozoa, Saprospiria, Methylacidiphilae,
    Calditrichae, Theionarchaea, Hydrogenophilalia, Udeonychophora,
    Priapulimorpha, Abditibacteria, Coprothermobacteria, Candidatus
    Fermentibacteria (class), Ichthyostraca, Palmophyllophyceae,
    Physodermatomycetes, Olpidiomycetes, Collemopsidiomycetes, Candidatus
    Sericytochromatia, Candidatus Riflemargulisbacteria, Endogonomycetes,
    Mucoromycetes, Umbelopsidomycetes, Mortierellomycetes, Kickxellomycetes,
    Dimargaritomycetes, Harpellomycetes, Asellariomycetes, Bartheletiomycetes,
    Zoopagomycetes
    '''
    classes = {i.strip() for i in classes_str.split(',')}
    # order exceptions that not end with "ales"
    order_exceptions = '''
    Labiatae, Gramineae, Kinetoplastida, Physariida, Haemosporida,
    Heterotrichida, Haptorida, Prorodontida, Haplosclerida, Actiniaria,
    Scleractinia, Pennatulacea, Semaeostomeae, Tricladida, Polycladida,
    Lecithoepitheliata, Strigeidida, Opisthorchiida, Cyclophyllidea,
    Heteronemertea, Monostilifera, Rhabditida, Ascaridida, Spirurida,
    Tylenchida, Trichinellida, Capitellida, Phyllodocida, Sabellida,
    Terebellida, Haplotaxida, Rhynchobdellida, Xenopneusta, Neogastropoda,
    Gymnosomata, Stylommatophora, Mytiloida, Arcoida, Ostreoida, Veneroida,
    Octopoda, Neoloricata, Anostraca, Decapoda, Euphausiacea, Amphipoda,
    Calanoida, Arguloida, Xiphosura, Scorpiones, Araneae, Ixodida, Odonata,
    Orthoptera, Phasmatodea, Coleoptera, Lepidoptera, Diptera, Hymenoptera,
    Mantodea, Siphonaptera, Neuroptera, Hemiptera, Porocephalida, Lingulida,
    Terebratulida, Comatulida, Velatida, Forcipulatida, Ophiurida, Cidaroida,
    Echinoida, Aspidochirotida, Molpadiida, Apodida, Dendrochirotida,
    Stolidobranchia, Petromyzontiformes, Myxiniformes, Rajiformes,
    Chimaeriformes, Lepidosireniformes, Coelacanthiformes, Acipenseriformes,
    Semionotiformes, Amiiformes, Elopiformes, Anguilliformes, Cypriniformes,
    Characiformes, Siluriformes, Gymnotiformes, Salmoniformes, Esociformes,
    Gadiformes, Batrachoidiformes, Lophiiformes, Atheriniformes, Perciformes,
    Pleuronectiformes, Beryciformes, Polypteriformes, Caudata, Anura,
    Gymnophiona, Testudines, Sphenodontia, Squamata, Casuariiformes,
    Rheiformes, Struthioniformes, Tinamiformes, Dinornithiformes,
    Apterygiformes, Anseriformes, Apodiformes, Caprimulgiformes,
    Charadriiformes, Ciconiiformes, Columbiformes, Coraciiformes,
    Cuculiformes, Falconiformes, Galliformes, Gruiformes, Passeriformes,
    Pelecaniformes, Phoenicopteriformes, Piciformes, Psittaciformes,
    Sphenisciformes, Turniciformes, Monotremata, Insectivora, Scandentia,
    Chiroptera, Primates, Cetacea, Sirenia, Proboscidea, Perissodactyla,
    Hyracoidea, Tubulidentata, Pholidota, Lagomorpha, Rodentia,
    Cheilostomatida, Dicyemida, Moniliformida, Arhynchobdellida, Cumacea,
    Mecoptera, Dermaptera, Mermithida, Plagiorchiida, Rhabdocoela, Cestida,
    Lobata, Poecilosclerida, Perkinsida, Choanoflagellida, Macroscelidea,
    Cyprinodontiformes, Pseudophyllidea, Gonorynchiformes, Rotaliida, Isopoda,
    Cephalobaenida, Archaeognatha, Ephemeroptera, Psocoptera, Strepsiptera,
    Thysanoptera, Trichoptera, Zygentoma, Zoraptera, Spirobolida,
    Spirostreptida, Salpida, Gaviiformes, Opisthocomiformes, Podicipediformes,
    Procellariiformes, Strigiformes, Carcharhiniformes, Hexanchiformes,
    Lamniformes, Orectolobiformes, Dermoptera, Tetraodontiformes, Zeiformes,
    Diadematoida, Arbacoida, Temnopleuroida, Clypeasteroida, Hymenostomatida,
    Clupeiformes, Nautilida, Dentaliida, Gadilida, Rhynchonellida,
    Aphragmophora, Phragmophora, Carnivora, Cephalodiscida, Rhabdopleurida,
    Bicosoecida, Peniculida, Loxodida, Armophorida, Euplotida, Stichotrichida,
    Vestibuliferida, Mesostigmata, Phymosomatoida, Cassiduloida,
    Eugregarinorida, Bivalvulida, Azygiida, Trichomonadida, Cyrtophorida,
    Ceriantharia, Stauromedusae, Rhizostomeae, Beroida, Protostomatida,
    Pygophora, Apygophora, Kentrogonida, Didelphimorphia, Paucituberculata,
    Microbiotheria, Dasyuromorphia, Diprotodontia, Notoryctemorphia,
    Peramelemorphia, Macrostomida, Proseriata, Nassulida, Polyopisthocotylea,
    Entodiniomorphida, Heterodontiformes, Alcyonacea, Schizopyrenida,
    Valvatida, Harpacticoida, Spatangoida, Chaetonotida, Solifugae,
    Scutigeromorpha, Scolopendromorpha, Lithobiomorpha, Echinothurioida,
    Julida, Osmeriformes, Osteoglossiformes, Echiuroinea, Mugiliformes,
    Solemyoida, Slopalinida, Siphonophorae, Trachymedusae, Narcomedusae,
    Schizomida, Opiliones, Acrasida, Synbranchiformes, Antipatharia,
    Corallimorpharia, Pterioida, Zoantharia, Polymorphida, Neoechinorhynchida,
    Clathrinida, Dorylaimida, Astrorhizida, Spionida, Nuculoida,
    Gyracanthocephala, Unionoida, Trigonioida, Stomiiformes, Cydippida,
    Raphidioptera, Megaloptera, Plecoptera, Embioptera, Diplogasterida,
    Galaxiiformes, Homalorhagida, Pseudoscorpiones, Ctenostomatida,
    Polyxenida, Proteocephalidea, Chaetodermatida, Limifossorida,
    Albuliformes, Vampyromorpha, Euryalida, Elasipodida, Spinulosida,
    Paxillosida, Isocrinida, Branchiobdellida, Musophagiformes, Trogoniformes,
    Trypanorhyncha, Tetraphyllidea, Echinorhynchida, Enoplida, Bucerotiformes,
    Coliiformes, Upupiformes, Ricinulei, Uropygi, Leptostraca,
    Grylloblattodea, Notostraca, Craterostigmomorpha, Geophilomorpha,
    Gigantorhynchida, Oligacanthorhynchida, Lyssacinosida, Acanthobdellida,
    Amblypygi, Glomerida, Pholadomyoida, Eunicida, Chaunacanthida,
    Symphyacanthida, Spumellaria, Holothyrida, Oxymonadida, Leptomyxida,
    Haplopharyngida, Bursovaginoidea, Thecideida, Aulopiformes,
    Myctophiformes, Multivalvulida, Monhysterida, Rhigonematida, Mononchida,
    Chromadorida, Oxyurida, Nudibranchia, Polydesmida, Siphonostomatoida,
    Desmodorida, Architaenioglossa, Amphionidacea, Bathynellacea, Anaspidacea,
    Mictacea, Tanaidacea, Mysida, Lophogastrida, Spelaeogriphacea,
    Thermosbaenacea, Eucoccidiorida, Beloniformes, Doliolida, Pyrosomata,
    Macrodasyida, Symphypleona, Neelipleona, Prolecithophora,
    Homosclerophorida, Percopsiformes, Polymixiiformes, Lampriformes,
    Trombidiformes, Sarcoptiformes, Cyclopoida, Platycopioida,
    Poecilostomatoida, Monstrilloida, Gelyelloida, Misophrioida,
    Mormonilloida, Podocopida, Myodocopida, Brachypoda, Diplostraca, Ploima,
    Arthracanthida, Lumbriculida, Agelasida, Phthiraptera, Blattodea,
    Leucosolenida, Coronatae, Helioporacea, Ptychodactiaria, Palpigradi,
    Bryophryida, Cyrtolophosidida, Colpodida, Bursariomorphida, Bryometopida,
    Bdellonemertea, Opilioacarida, Telestacea, Clevelandellida, Plagiotomida,
    Ophidiiformes, Cyclostomatida, Hexactinosida, Limoida, Myzostomida,
    Adinetida, Philodinida, Flosculariacea, Pedinoida, Brisingida, Pectinoida,
    Amphilinidea, Spathebothriidea, Caryophyllidea, Haplobothriidea,
    Diphyllidea, Lecanicephalidea, Tetrabothriidea, Rhinebothriidea,
    Araeolaimida, Craniida, Pedunculata, Sessilia, Dendrogastrida, Laurida,
    Mystacocaridida, Nectiopoda, Platycopida, Nippotaeniidea, Myliobatiformes,
    Pristiformes/Rhiniformes group, Torpediniformes, Echinorhiniformes,
    Ceratodontiformes, Galbuliformes, Callipodida, Chordeumatida,
    Sphaerotheriida, Polyzoniida, Playtdesmida, Dendroceratida,
    Dictyoceratida, Filospermoidea, Gyrocotylidea, Limnomedusae,
    Syngnathiformes, Litobothriidea, Peripodida, Argentiniformes, Stemonitida,
    Platyctenida, Thalassocalycida, Ateleopodiformes, Notacanthiformes,
    Parachela, Apochela, Thecosomata, Exogenida, Evaginogenida, Endogenida,
    Sorogenida, Euglyphida, Halocyprida, Enterogona, Cercomonadida,
    Thaumatomonadida, Chordodea, Gordea, Diplonemida, Mantophasmatodea,
    Pleurostomatida, Cyclotrichida, Philasterida, Pleuronematida,
    Thigmotrichida, Dermocystida, Ichthyophonida, Kiitrichida, Tintinnida,
    Choreotrichida, Trichiida, Triplonchida, Isolaimida, Licnophorida,
    Aspidosiphonidormes, Phascolosomatiformes, Littorinimorpha,
    Protoheterotrichida, Bourgueticrinida, Cyrtocrinida, Exogemmida,
    Neogregarinorida, Phaeosphaerida, Phaeodendrida, Phaeocystida,
    Chlamydodontida, Dysteriida, Glomeridesmida, Stemmiulida, Acochlidiacea,
    Amphidiscosida, Potamodrilidae, Myoida, Astomatida, Nassellaria,
    Flabelligerida, Squaliformes, Murrayonida, Lithonida, Trichonymphida,
    Spirotrichonymphida, Cyclorhagida, Himatismenida, Arcellinida,
    Akentrogonida, Cathetocephalidea, Actinulida, Holasteroida, Phaeogromida,
    Phaeoconchida, Philodinavidae, Synhymeniida, Notomyotida, Tubuliporida,
    Monophragmophora, Pantopoda, Desmoscolecida, Pseudoamphisiellida,
    Anthoathecata, Odontostomatida, Siphonophorida, Agamococcidiorida,
    Grossglockneriida, Urostylida, Leptothecata, Microthoracida,
    Malacovalvulida, Apostomatida, Nanaloricida, Sepiida, Sepiolida,
    Spirulida, Teuthida, Telonemida, Cryomonadida, Polystilifera, Carybdeida,
    Chirodropida, Baerida, Pleurobranchomorpha, Moniligastrida,
    Sporadotrichida, Hyocrinida, Arthrotardigrada, Echiniscoidea, Opheliida,
    Entomobryomorpha, Poduromorpha, Echinothuroida, Millericrinida,
    Tritrichomonadida, Honigbergiellida, Cristamonadida, Hypotrichomonadida,
    Golfingiida, Nuculanoida, Crustaceacida, Pilosa, Cingulata, Cycloneritida,
    Lituolida, Spirillinida, Sticholonchida, Neomeniamorpha, Pholidoskepia,
    Cavibelonia, Bacteroidetes Order II. Incertae sedis,Bacteroidetes Order
    IV. Incertae sedis, Holacanthida, Muspiceida, Rhynchodida, Hypocomatida,
    Rapaza, Miliolida, Diphyllobothriidea, Rigifilida, Glissomonadida,
    Plagiopylida, Archigregarinorida, Crocodylia, Fecampiida, Picomonadida,
    Bothriocephalidea, Vampyrellida, Carterinida, Phyllobothriidea,
    Onchoproteocephalidea, Salenioida, Collodaria, Longamoebia,
    Hiodontiformes, Alepocephaliformes, Lepidogalaxiiformes, Stylephoriformes,
    Kurtiformes, Gobiiformes, Scombriformes, Anabantiformes, Istiophoriformes,
    Carangiformes, Cichliformes, Pholidichthyiformes, Blenniiformes,
    Uranoscopiformes, Labriformes, Lobotiformes, Ephippiformes, Spariformes,
    Acanthuriformes, Pempheriformes, Centrarchiformes, Holocentriformes,
    Chaetodontiformes, Notoungulata, Litopterna, Lucinoida, Phaeocalpida,
    Phaeogymnocellida, Collothecaceae, Protura, Robertinida, Prostomatida,
    Chondrillida, Verongiida, Axinellida, Biemnida, Bubarida, Clionaida,
    Desmacellida, Merliida, Polymastiida, Sphaerocladina, Spongillida,
    Suberitida, Tethyida, Tetractinellida, Trachycladida, Aulocalycoida,
    Cariamiformes, Enchytraeida, Sessilida, Mobilida, Caproiformes,
    Priacanthiformes, Gerreiformes, Lutjaniformes, Dioctophymatida,
    Benthimermithida, Plectida, Euonychophora, Priapulimorphida,
    Tetramerocerata, Rhaptothyreida, Sterrofustia, Aquavolonida, Micropygoida,
    Metopida, Hirudinida
    '''
    orders = {i.strip() for i in order_exceptions.split(',')}

    def get_taxon(taxon_str):
        """
        Get taxon info based on suffix and list from NCBI taxonomy database.
        """
        # kingdom|phylum|class|order|family|organims(genus|species)
        # add my_ prefix to avoid conflict of "class"
        my_kingdom = ''
        my_phylum = ''
        my_class = ''
        my_order = ''
        my_family = ''
        for item in taxon_str:
            if item in superkingdoms:
                my_kingdom = item
            # mix superkingdom and kingdom to reduce name length
            elif item in kingdoms:
                my_kingdom = item
            elif item in phyla:
                my_phylum = item
            elif item in classes:
                my_class = item
            if item.endswith('ales') or item in orders:
                my_order = item
            elif item.endswith('aceae') or item.endswith('idae'):
                my_family = item
        # get fake class for plant
        if my_phylum == 'Streptophyta' and my_class == '':
            last_phyta = ''
            for i in taxon_str:
                if i.endswith('phyta'):
                    last_phyta = i
            try:
                my_class = taxon_str[taxon_str.index(last_phyta) + 1]
            except IndexError:
                my_class = ''
        return (my_kingdom, my_phylum, my_class, my_order, my_family)

    # put raw fasta into root of output folder, so not to use clean_path
    raw_fasta = join_path(arg.out, splitext(basename(gbfile))[0] + '.fasta')
    handle_raw = open(raw_fasta, 'w', encoding='utf-8')
    wrote_by_gene = set()
    wrote_by_name = set()
    log.info('Divide {}.'.format(gbfile))
    sequence_id = 0
    for record in clean_gb(gbfile):
        # only accept gene, product, and spacer in misc_features.note
        taxon_str = record.annotations['taxonomy']
        kingdom, phylum, class_, order, family = get_taxon(taxon_str)
        organism = record.annotations['organism'].replace(' ', '_')
        genus, *species = organism.split('_')
        # species name may contain other characters
        taxon = '{}|{}|{}|{}|{}|{}|{}'.format(kingdom, phylum, class_,
                                              order, family, genus,
                                              '_'.join(species))
        try:
            sequence_id = accession = record.annotations['accessions'][0]
        except KeyError:
            sequence_id = accession = ''
        try:
            specimen = record.features[0].qualifiers[
                'specimen_voucher'][0].replace(' ', '_')
        except (IndexError, KeyError):
            specimen = ''
        whole_seq = record.seq
        feature_name = []
        genes = []

        for feature in record.features:
            if feature.type == 'source':
                continue
            name, feature_type = get_feature_name(feature, arg)
            # skip unsupport feature
            if name is None:
                continue
            if len(name) > arg.max_name_len:
                log.warning('Too long name: {}. Truncated.'.format(name))
                name = name[:arg.max_name_len] + '...'
            # skip abnormal annotation
            if len(feature) > arg.max_seq_len:
                log.warning('Skip abnormal annotaion of {}'
                            '(Accession {}).'.format(name, accession))
                continue
            if feature_type == 'gene':
                genes.append([name, feature])
            feature_name.append(name)
            sequence_id = '>' + '|'.join([name, taxon, accession,
                                          specimen])
            wrote = write_seq(name, sequence_id, feature, whole_seq,
                              arg.by_gene_folder, arg)
            wrote_by_gene.add(wrote)

        # extract spacer
        spacers = get_spacer(genes)
        for spacer in spacers:
            if len(spacer) > arg.max_seq_len:
                log.warning('Spacer {} too long (Accession {}).'.format(
                    spacer.id, accession))
                continue
            sequence_id = '>' + '|'.join([spacer.id, taxon,
                                          accession, specimen])
            wrote = write_seq(spacer.id, sequence_id, spacer, whole_seq,
                              arg.by_gene_folder, arg)
            wrote_by_gene.add(wrote)
        # write to group_by name, i.e., one gb record one fasta
        if 'ITS' in feature_name:
            name_str = 'ITS'
        elif len(feature_name) >= 4:
            name_str = '{}-...-{}'.format(feature_name[0], feature_name[-1])
        elif len(feature_name) == 0:
            name_str = 'Unknown'
        else:
            name_str = '-'.join(feature_name)
        # directly use genome type as name
        if arg.organelle is not None:
            name_str = '{}_genome'.format(arg.organelle)
        record.id = '|'.join([name_str, taxon, accession, specimen])
        record.description = ''
        filename = join_path(arg.by_name_folder, name_str + '.fasta')
        with open(filename, 'a', encoding='utf-8') as out:
            SeqIO.write(record, out, 'fasta')
            wrote_by_name.add(filename)
        # write raw fasta
        SeqIO.write(record, handle_raw, 'fasta')

    # skip analyze of Unknown.fasta
    unknown = join_path(arg.by_name_folder, 'Unknown.fasta')
    if unknown in wrote_by_name:
        log.info('Skip Unknown.fasta')
        wrote_by_name.remove(unknown)
    log.info('Divide finished.')
    return list(wrote_by_gene), list(wrote_by_name)


def uniq(files, arg):
    """
    Remove redundant sequences of same species.
    """
    uniq_files = []
    for fasta in files:
        info = defaultdict(lambda: list())
        keep = dict()
        index = 0
        for record in SeqIO.parse(fasta, 'fasta'):
            # gene|kingdom|phylum|class|order|family|genus|species|specimen
            if '|' in record.id:
                name = ' '.join(record.id.split('|')[6:8])
            else:
                name = record.id
            length = len(record)
            # skip empty file
            if length != 0:
                info[name].append([index, length])
            index += 1
        if arg.uniq == 'first':
            # keep only the first record
            keep = {info[i][0][0] for i in info}
        elif arg.uniq == 'longest':
            for i in info:
                info[i] = sorted(info[i], key=lambda x: x[1], reverse=True)
            keep = {info[i][0][0] for i in info}
        elif arg.uniq == 'random':
            for i in info:
                info[i] = choice(info[i])
            keep = {info[i][0] for i in info}
        new = clean_path(fasta, arg) + '.uniq'
        with open(new, 'w', encoding='utf-8') as out:
            for idx, record in enumerate(SeqIO.parse(fasta, 'fasta')):
                if idx in keep:
                    SeqIO.write(record, out, 'fasta')
        log.info('{} of {} records were kept in {}'.format(len(keep), index,
                                                           fasta))
        uniq_files.append(new)
    return uniq_files


def align(files, arg):
    """
    Calls mafft to align sequences.
    """
    log.info('Align sequences.')
    result = []
    # get available CPU cores
    cores = max(1, cpu_count() - 1)
    for fasta in files:
        log.info('Aligning {}.'.format(fasta))
        out = clean_path(fasta, arg) + '.aln'
        with open(devnull, 'w', encoding='utf-8') as f:
            # if computer is good enough, "--genafpair" is recommended
            _ = ('mafft --thread {} --reorder --quiet --adjustdirection '
                 '{} > {}'.format(cores, fasta, out))
            m = run(_, shell=True, stdout=f, stderr=f)
        if m.returncode == 0:
            result.append(out)
        else:
            # ignore empty result
            pass
    log.info('Alignment finished.')
    for i in glob('_order*'):
        remove(i)
    log.info('{} of {} files were successfully aligned.'.format(len(result),
                                                                len(files)))
    return result


def prepare(aln_fasta, arg):
    """
    Given fasta format alignment filename, return a numpy array for sequence:
    Generate fasta file without gap for makeblastdb, return file name.
    Faster and use smaller mem :)
    """
    data = []
    record = ['id', 'sequence']
    no_gap = StringIO()
    with open(aln_fasta, 'r', encoding='utf-8') as raw:
        for line in raw:
            no_gap.write(line.replace('-', ''))
            if line.startswith('>'):
                data.append([record[0], ''.join(record[1:])])
                # remove ">" and CRLF
                name = line.strip('>\r\n')
                record = [name, '']
            else:
                record.append(line.strip().upper())
        # add last sequence
        data.append([record[0], ''.join(record[1:])])
    # skip head['id', 'seq']
    data = data[1:]
    # check sequence length
    length_check = [len(i[1]) for i in data]
    if len(set(length_check)) != 1:
        log.error('{} does not have uniform width!'.format(aln_fasta))
        return None, None, None

    # Convert List to numpy array.
    # order 'F' is a bit faster than 'C'
    # new = np.hstack((name, seq)) -> is slower
    name = np.array([[i[0]] for i in data], dtype=np.bytes_)
    sequence = np.array([list(i[1]) for i in data], dtype=np.bytes_, order='F')

    if name is None:
        log.error('Bad fasta file {}.'.format(aln_fasta))
        name = None
    # try to avoid makeblastdb error
    no_gap.seek(0)
    SeqIO.convert(no_gap, 'fasta', arg.db_file, 'fasta')
    no_gap.close()
    return name, sequence, arg.db_file


def count_base(alignment, rows, columns):
    """
    Given alignment numpy array, count cumulative frequency of base in each
    column (consider ambiguous base and "N", "-" and "?", otherwise omit).
    Return [[float, float, float, float, float, float, float]] for
    [A, T, C, G, N, GAP, OTHER].
    """
    frequency = []
    for index in range(columns):
        base, counts = np.unique(alignment[:, [index]], return_counts=True)
        count_dict = {b'A': 0, b'C': 0, b'G': 0, b'T': 0, b'M': 0, b'R': 0,
                      b'W': 0, b'S': 0, b'Y': 0, b'K': 0, b'V': 0, b'H': 0,
                      b'D': 0, b'B': 0, b'X': 0, b'N': 0, b'-': 0, b'?': 0}
        count_dict.update(dict(zip(base, counts)))
        a = (count_dict[b'A'] +
             (count_dict[b'D'] + count_dict[b'H'] + count_dict[b'V']) / 3 +
             (count_dict[b'M'] + count_dict[b'R'] + count_dict[b'W']) / 2)
        t = (count_dict[b'T'] +
             (count_dict[b'B'] + count_dict[b'H'] + count_dict[b'D']) / 3 +
             (count_dict[b'K'] + count_dict[b'W'] + count_dict[b'Y']) / 2)
        c = (count_dict[b'C'] +
             (count_dict[b'B'] + count_dict[b'H'] + count_dict[b'V']) / 3 +
             (count_dict[b'M'] + count_dict[b'S'] + count_dict[b'Y']) / 2)
        g = (count_dict[b'G'] +
             (count_dict[b'B'] + count_dict[b'D'] + count_dict[b'V']) / 3 +
             (count_dict[b'K'] + count_dict[b'R'] + count_dict[b'S']) / 2)
        gap = count_dict[b'-']
        n = count_dict[b'N'] + count_dict[b'X'] + count_dict[b'?']
        other = rows - a - t - c - g - gap - n
        frequency.append([a, t, c, g, n, gap, other])
    return frequency


def get_quality(data, rows):
    """
    Calculate quality score.
    """
    # use fastq-illumina format
    max_q = 62
    factor = max_q / rows
    # use min to avoid KeyError
    quality_value = [min(max_q, int(i * factor)) - 1 for i in data]
    return quality_value


def get_resolution(alignment, start, end, fast=False):
    """
    Given alignment (2d numpy array), location of fragment(start and end, int,
    start from zero, exclude end),
    return gap ratio, resolution, entropy, Pi, tree value and average terminal
    branch length.
    """
    subalign = alignment[:, start:end]
    rows, columns = subalign.shape
    max_h = np.log2(rows)
    total = rows * columns
    gap_ratio = 0
    resolution = 0
    entropy = 0
    pi = 0
    tree_value = 0
    avg_terminal_branch_len = 0
    # index error
    if columns == 0:
        return (gap_ratio, resolution, entropy, pi, tree_value,
                avg_terminal_branch_len)
    gap_ratio = len(subalign[subalign == b'-']) / total
    item, count = np.unique(subalign, return_counts=True, axis=0)
    resolution = len(count) / rows
    # normalized entropy
    entropy = 0
    for j in count:
        p_j = j / rows
        log2_p_j = np.log2(p_j)
        entropy += log2_p_j * p_j
    entropy = -1 * entropy / max_h
    # Nucleotide diversity (pi)
    m = columns
    n = rows
    sum_d_ij = 0
    for i in range(n):
        d_ij = np.sum(subalign[i] != subalign[(i + 1):])
        sum_d_ij += d_ij
    pi = (2 / (n * (n - 1)) * sum_d_ij) / m
    # tree value
    aln_file = '{}-{}.aln.tmp'.format(start, end)

    def clean():
        for _ in glob(aln_file + '*'):
            remove(_)
    if not fast:
        if rows < 4:
            return (gap_ratio, resolution, entropy, pi, tree_value,
                    avg_terminal_branch_len)
        with open(aln_file, 'wb') as aln:
            for index, row in enumerate(alignment[:, start:end]):
                aln.write(b'>' + str(index).encode('utf-8') + b'\n' + b''.join(
                    row) + b'\n')
        with open(devnull, 'w', encoding='utf-8') as f:
            iqtree = run('iqtree -s {} -m JC -fast -czb'.format(aln_file),
                         stdout=f, stderr=f, shell=True)
        # just return 0 if there is error
        if iqtree.returncode != 0:
            log.info('Too much gap in the region {}-{} bp.'.format(start, end))
            clean()
        else:
            tree = Phylo.read(aln_file + '.treefile', 'newick')
            # skip the first empty node
            internals = tree.get_nonterminals()[1:]
            terminals = tree.get_terminals()
            sum_terminal_branch_len = sum([i.branch_length for i in terminals])
            # miss stdev, to be continued
            avg_terminal_branch_len = sum_terminal_branch_len / len(terminals)
            tree_value = len(internals) / len(terminals)
            clean()
    return (gap_ratio, resolution, entropy, pi, tree_value,
            avg_terminal_branch_len)


def generate_consensus(base_cumulative_frequency, coverage_percent,
                       rows, output):
    """
    Given count info of bases, return consensus(PrimerWithInfo).
    """
    def get_ambiguous_dict():
        data = dict(zip(ambiguous_data.values(), ambiguous_data.keys()))
        # 2:{'AC': 'M',}
        data_with_len = defaultdict(lambda: dict())
        for i in data:
            data_with_len[len(i)][i] = data[i]
        return data_with_len

    ambiguous_dict = get_ambiguous_dict()
    most = []
    coverage = rows * coverage_percent

    limit = coverage / len('ATCG')
    for location, column in enumerate(base_cumulative_frequency):
        finish = False
        # "*" for others
        value = dict(zip(list('ATCGN-*'), column))

        base = 'N'
        if value['N'] >= limit:
            count = value['N']
            most.append([location, base, count])
            continue
        sum_gap = value['-'] + value['*']
        if sum_gap >= limit:
            base = '-'
            count = sum_gap
            most.append([location, base, count])
            continue
        # 1 2 3 4
        for length in ambiguous_dict:
            # A T CG CT ACG CTG ATCG
            for key in ambiguous_dict[length]:
                count = 0
                for letter in list(key):
                    if finish:
                        break
                    count += value[letter]
                    if count >= coverage:
                        base = ambiguous_dict[length][key]
                        finish = True
                        most.append([location, base, count])
    quality_raw = [i[2] for i in most]
    consensus = PrimerWithInfo(start=1, seq=''.join([i[1] for i in most]),
                               quality=get_quality(quality_raw, rows))
    SeqIO.write(consensus, output, 'fastq')
    return consensus


def get_good_region(index, seq_count, arg):
    """
    Return regions marked for finding primers. Because of alignment gap, PCR
    product may smaller than given length limitation.
    """
    n = arg.max_product - arg.min_product
    good_region = set()
    for i, j in zip(index, seq_count):
        if j >= arg.resolution:
            good_region.update(range(i - arg.max_primer, i - n))
            good_region.update(range(i + arg.min_product,
                                     i - arg.max_primer + arg.max_product))
    return good_region


def find_continuous(consensus, good_region, min_len):
    """
    Given PrimerWithInfo, good_region, min_len
    Return consensus with features.
    """
    skip = ('N', '-')
    start = 0
    for index, base in enumerate(consensus.sequence[:-min_len]):
        if base in skip or index not in good_region:
            if (index - start) >= min_len:
                consensus.features.append(SeqFeature(FeatureLocation(
                    start, index), type='continuous', strand=1))
            start = index + 1
    return consensus


def find_primer(consensus, arg):
    """
    Find suitable primer in given consensus with features labeled as candidate
    primer, return list of PrimerWithInfo, consensus.
    """
    # repeat no more than 5 times
    poly = re.compile(r'([ATCG])\1\1\1\1')
    tandem = re.compile(r'([ATCG]{2})\1\1\1\1')

    def is_good_primer(primer):
        # use re and primer3 to check weather it's good primer
        # ref1. http://www.premierbiosoft.com/tech_notes/PCR_Primer_Design.html
        ambiguous_base = len(primer)
        for i in list('ATCG'):
            ambiguous_base -= primer.seq.count(i)
        if ambiguous_base > max_ambiguous:
            return False
        if re.search(poly, str(primer.seq)) is not None:
            primer.detail = 'Poly(NNNNN) structure found'
            return False
        if re.search(tandem, str(primer.seq)) is not None:
            primer.detail = 'Tandom(NN*5) exist'
            return False
        primer.hairpin_tm = calc_ambiguous_seq(calcHairpinTm, primer.seq)
        primer.tm = primer.annotations['tm'] = calc_ambiguous_seq(calcTm,
                                                                  primer.seq)
        # primer3.calcHairpin or calcHomodimer usually return structure found
        # with low Tm. Here we compare structure_tm with sequence tm
        if primer.hairpin_tm >= primer.tm:
            primer.detail = 'Hairpin found'
            return False
        primer.homodimer_tm = calc_ambiguous_seq(calcHomodimerTm, primer.seq)
        if primer.homodimer_tm >= primer.tm:
            primer.detail = 'Homodimer found'
            return False
        return True

    primers = []
    min_len = arg.min_primer
    max_len = arg.max_primer
    max_ambiguous = arg.ambiguous_base_n
    # skip good_region
    continuous = consensus.features
    for feature in continuous:
        fragment = feature.extract(consensus)
        len_fragment = len(fragment)
        for begin in range(len_fragment - max_len):
            for p_len in range(min_len, max_len + 1):
                start = feature.location.start + begin
                primer = consensus[start:start + p_len]
                if is_good_primer(primer):
                    consensus.features.append(SeqFeature(
                        FeatureLocation(start, start + p_len),
                        type='primer', strand=1))
                    primer.start = start
                    primer.update_id()
                    primers.append(primer)
    return primers, consensus


def count_and_draw(alignment, arg):
    """
    Given alignment(numpy array), calculate Shannon Index based on
    www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/shannon.htm
    Return lists of observed resolution, shannon index, Pi, tree resolution,
    average terminal branch length and index.
    Draw sliding-window figure.
    All calculation excludes primer sequence.
    """
    output = join_path(arg.out, basename(arg.out_file).split('.')[0])
    rows, columns = alignment.shape
    min_primer = arg.min_primer
    max_product = arg.max_product
    step = arg.step
    if rows < 4:
        log.warning('Less than 3 sequence. Tree inference will be skipped.')
    # gap_list, r_list, h_list, pi_list, t_list : count, normalized entropy,
    # Pi and tree value
    gap_ratio_list = []
    entropy_list = []
    avg_branch_len_list = []
    pi_list = []
    observed_res_list = []
    tree_res_list = []
    index = []
    max_plus = max_product - min_primer * 2
    max_range = columns - max_product
    handle = open(output + '.variance.tsv', 'w', encoding='utf-8')
    # iqtree blmin is 1e-6
    fmt = '{},{:.2%},{:.6f},{:.6f},{:.6f},{:.6f},{:.6f}\n'
    handle.write('Index,GapRatio,Resolution,Entropy,Pi,'
                 'TreeValue,AvgTerminalBranchLen\n')
    for i in range(0, max_range, step):
        # exclude primer sequence
        values = get_resolution(alignment, i, i + max_plus, arg.fast)
        handle.write(fmt.format(i, *values))
        gap_ratio, resolution, entropy, pi, tree_value, avg_branch_len = values
        gap_ratio_list.append(gap_ratio)
        observed_res_list.append(resolution)
        entropy_list.append(entropy)
        pi_list.append(pi)
        tree_res_list.append(tree_value)
        avg_branch_len_list.append(avg_branch_len)
        index.append(i)

    plt.style.use('seaborn-colorblind')
    # how to find optimized size?
    fig, ax1 = plt.subplots(figsize=(15 + len(index) // 5000, 10))
    plt.title('Variance of {} (sample={}, window={} bp, step={} bp)\n'.format(
        basename(arg.out_file).split('.')[0], rows, max_product, step))
    plt.xlabel('Base')
    ax1.yaxis.set_ticks(np.linspace(0, 1, num=11))
    ax1.set_ylabel('Resolution & Shannon Equitability Index')
    ax1.plot(index, entropy_list, label='Shannon Equitability Index',
             alpha=0.8)
    ax1.plot(index, observed_res_list, label='Observed Resolution', alpha=0.8)
    ax1.plot(index, gap_ratio_list, label='Gap Ratio', alpha=0.8)
    # different ytick
    ax2 = ax1.twinx()
    ax2.set_ylabel(r'$\pi$', rotation=-90, labelpad=20)
    # plt.xticks(np.linspace(0, max_range, 21))
    if not arg.fast:
        ax1.plot(index, tree_res_list, label='Tree Resolution', alpha=0.8)
        ax2.plot(index, avg_branch_len_list, linestyle='--',
                 label='Average Terminal Branch Length', alpha=0.8)
        ax2.set_ylabel(r'$\pi$ and Average Branch Length', rotation=-90,
                       labelpad=20)
    ax1.legend(loc='lower left')
    ax2.plot(index, pi_list, 'k--', label=r'$\pi$', alpha=0.8)
    ax2.legend(loc='upper right')
    plt.savefig(output + '.pdf')
    plt.savefig(output + '.png')
    plt.close()
    handle.close()
    return observed_res_list, index


def parse_blast_tab(filename):
    """
    Parse BLAST result (tab format).
    """
    query = []
    with open(filename, 'r', encoding='utf-8') as raw:
        for line in raw:
            if line.startswith('# BLAST'):
                yield query
                query.clear()
            elif line.startswith('#'):
                pass
            else:
                query.append(BlastResult(line))


def validate(primer_candidate, db_file, n_seqs, arg):
    """
    Do BLAST. Parse BLAST result. Return list of PrimerWithInfo which passed
    the validation.
    """
    EVALUE = 1e-2
    query_file = arg.out_file + '.candidate.fasta'
    query_file_fastq = arg.out_file + '.candidate.fastq'
    # SeqIO.write fasta file directly is prohibited. have to write fastq at
    # first.
    with open(query_file_fastq, 'w', encoding='utf-8') as _:
        SeqIO.write(primer_candidate, _, 'fastq')
    SeqIO.convert(query_file_fastq, 'fastq', query_file, 'fasta')
    # build blast db
    with open(devnull, 'w', encoding='utf-8') as f:
        _ = run('makeblastdb -in {} -dbtype nucl'.format(db_file),
                shell=True, stdout=f)
        if _.returncode != 0:
            log.critical('Failed to run makeblastdb. Skip BLAST.')
            return []
    # BLAST
    blast_result_file = 'blast.result.tsv'
    fmt = 'qseqid sseqid qseq nident mismatch score qstart qend sstart send'
    cmd = Blast(num_threads=max(1, cpu_count() - 1),
                query=query_file,
                db=db_file,
                task='blastn-short',
                evalue=EVALUE,
                max_hsps=1,
                outfmt='"7 {}"'.format(fmt),
                out=blast_result_file)
    # hide output
    cmd()
    log.info('BLAST finished.')
    log.info('Parsing BLAST result.')
    blast_result = dict()
    # because SearchIO.parse is slow, use parse_blast_result()
    for query in parse_blast_tab(blast_result_file):
        if len(query) == 0:
            continue
        sum_bitscore_raw = 0
        sum_mismatch = 0
        good_hits = 0
        mid_loc = dict()
        hit = query[0]
        for hit in query:
            min_positive = len(hit.query_seq) - arg.mismatch
            hsp_bitscore_raw = hit.bitscore_raw
            positive = hit.ident_num
            mismatch = hit.mismatch_num
            loc = average([hit.hit_start, hit.hit_end])
            if positive >= min_positive and mismatch <= arg.mismatch:
                sum_bitscore_raw += hsp_bitscore_raw
                sum_mismatch += mismatch
                good_hits += 1
                # middle location of primer, the difference of two mid_loc
                # approximately equals to the length of amplified fragment.
                mid_loc[hit.hit_id] = loc
        coverage = good_hits / n_seqs
        if coverage >= arg.coverage:
            blast_result[hit.query_id] = {
                'coverage': coverage,
                'avg_bitscore': sum_bitscore_raw / good_hits,
                'avg_mismatch': sum_mismatch / good_hits,
                'mid_loc': mid_loc}
    primer_verified = []
    for primer in primer_candidate:
        i = primer.id
        if i in blast_result:
            primer.coverage = blast_result[i]['coverage']
            primer.avg_bitscore = blast_result[i]['avg_bitscore']
            primer.mid_loc = blast_result[i]['mid_loc']
            primer.avg_mismatch = blast_result[i]['avg_mismatch']
            primer.update_id()
            primer_verified.append(primer)
    primer_verified.sort(key=lambda x: x.start)
    # clean makeblastdb files
    for i in glob(db_file + '*'):
        remove(i)
    remove(blast_result_file)
    log.info('Parse finished.')
    return primer_verified


def pick_pair(primers, alignment, arg):
    """
    Pick primer pairs passed the validation and its product length fulfill the
    requirement.
    """
    pairs = []
    for n_left, left in enumerate(primers):
        # convert mid_loc to 5' location
        # use int to speedup, comparing of float seems slow
        location = int(left.avg_mid_loc - len(left) / 2)
        begin = location + arg.min_product
        # fragment plus one primer = max_product length
        end = location + arg.max_product - len(left)
        # create [] is faster than list()
        cluster = []
        for right in primers[(n_left + 1):]:
            if right.avg_mid_loc < begin:
                continue
            if right.avg_mid_loc > end:
                break
            pair = Pair(left, right, alignment)
            if pair.coverage < arg.coverage:
                continue
            cluster.append(pair)
        cluster.sort(key=lambda x: x.score, reverse=True)
        # only keep top n for each primer cluster
        pairs.extend(cluster[:arg.top_n])
    if len(pairs) == 0:
        return []
    # remove close located primers
    less_pairs = []
    cluster = [pairs[0], ]
    pairs.sort(key=lambda x: x.start)
    for index in range(1, len(pairs)):
        if pairs[index].start - pairs[index - 1].start < arg.min_primer:
            cluster.append(pairs[index])
        else:
            cluster.sort(key=lambda x: x.score, reverse=True)
            less_pairs.extend(cluster[:arg.top_n])
            cluster.clear()
    cluster.sort(key=lambda x: x.score, reverse=True)
    less_pairs.extend(cluster[:arg.top_n])
    log.info('{} pairs of redundant primers were removed.'.format(
        len(pairs) - len(less_pairs)))
    good_pairs = []
    for i in less_pairs:
        i.add_info(alignment)
        if i.resolution >= arg.resolution:
            good_pairs.append(i)
    good_pairs.sort(key=lambda x: x.score, reverse=True)
    log.info('Successfully found validated primers.')
    return good_pairs[:arg.top_n]


def analyze(fasta, arg):
    """
    Primer design pipeline.
    Return bool for success or not.
    """
    # read from fasta, generate new fasta for makeblastdb
    name, alignment, db_file = prepare(fasta, arg)
    if name is None:
        log.info('Invalid fasta file {}.'.format(fasta))
        return False
    rows, columns = alignment.shape
    # generate consensus
    base_cumulative_frequency = count_base(alignment, rows, columns)
    log.info('Generate consensus of {}.'.format(fasta))
    consensus = generate_consensus(base_cumulative_frequency, arg.coverage,
                                   rows, arg.out_file + '.consensus.fastq')
    log.info('Evaluate whole alignment of {}.'.format(fasta))
    # a_ : alignment
    (a_gap_ratio, a_observed_res, a_entropy, a_pi, a_tree_res,
     a_branch_len) = get_resolution(alignment, 0, columns)
    log.info('Gap ratio: {}'.format(a_gap_ratio))
    log.info('Observed resolution: {}'.format(a_observed_res))
    log.info('Normalized Shannon Index: {}'.format(a_entropy))
    log.info('Pi: {}'.format(a_pi))
    log.info('Tree resolution: {}'.format(a_tree_res))
    log.info('Average terminal branch length: {}'.format(a_branch_len))
    summary = join_path(arg.out, 'Loci.csv')
    if not exists(summary):
        with open(summary, 'w', encoding='utf-8') as s:
            s.write('Loci,Samples,Length,GapRatio,ObservedResolution,'
                    'TreeResolution,ShannonIndex,AvgTerminalBranchLen,Pi\n')
            s.write('{},{},{},{:.2%},{:.6f},{:.6f},{:.6f},{:.6f},{:.6f}'
                    '\n'.format(basename(fasta), rows, columns, a_gap_ratio,
                                a_observed_res, a_tree_res, a_entropy,
                                a_branch_len, a_pi))
    else:
        with open(summary, 'a', encoding='utf-8') as s:
            s.write('{},{},{},{:.2%},{:.6f},{:.6f},{:.6f},{:.6f},{:.6f}'
                    '\n'.format(basename(fasta), rows, columns, a_gap_ratio,
                                a_observed_res, a_tree_res, a_entropy,
                                a_branch_len, a_pi))
    # exit if resolution lower than given threshold.
    if a_observed_res < arg.resolution:
        log.warning('Observed resolution is too low.')
        return False
    log.info('Start the sliding-window scan.')
    observed_res_list, index = count_and_draw(alignment, arg)
    log.info('Evaluation finished.')
    # stop if do not want to design primer
    if arg.stop == 2:
        return True
    log.info('Start finding primers.')
    log.info('Mark region for finding primer.')
    good_region = get_good_region(index, observed_res_list, arg)
    log.info('Finding candidate primers.')
    consensus = find_continuous(consensus, good_region, arg.min_primer)
    primer_candidate, consensus = find_primer(consensus, arg)
    if len(primer_candidate) == 0:
        log.warning('Cannot find primer candidates. '
                    'Please consider to loose options.')
        return True
    log.info('Found {} candidate primers.'.format(len(primer_candidate)))
    log.info('Validate with BLAST. May be slow.')
    primer_verified = validate(primer_candidate, db_file, rows, arg)
    if len(primer_verified) == 0:
        log.warning('All candidates failed on validation. '
                    'Please consider to loose options.')
        return True
    log.info('Picking primer pairs.')
    pairs = pick_pair(primer_verified, alignment, arg)
    if len(pairs) == 0:
        log.warning('Cannot find suitable primer pairs. '
                    'Please consider to loose options.')
        return True
    log.info('Output the result.')
    locus = basename(arg.out_file).split('.')[0]
    csv_title = ('Locus,Score,Samples,AvgProductLength,StdEV,'
                 'MinProductLength,MaxProductLength,'
                 'Coverage,Resolution,TreeValue,AvgTerminalBranchLen,Entropy,'
                 'LeftSeq,LeftTm,LeftAvgBitscore,LeftAvgMismatch,'
                 'RightSeq,RightTm,RightAvgBitscore,RightAvgMismatch,'
                 'DeltaTm,AlnStart,AlnEnd,AvgSeqStart,AvgSeqEnd\n')
    style = ('{},{:.2f},{},{:.0f},{:.0f},{},{},'
             '{:.2%},{:.2%},{:.6f},{:.6f},{:.6f},'
             '{},{:.2f},{:.2f},{:.2f},'
             '{},{:.2f},{:.2f},{:.2f},'
             '{:.2f},{},{},{},{}\n')
    out1 = open(join_path(arg.out, locus) + '.primer.fastq', 'w',
                encoding='utf-8')
    out2 = open(join_path(arg.out, locus) + '.primer.csv', 'w',
                encoding='utf-8')
    # write primers to one file
    out3_file = join_path(arg.out, 'Primers.csv')
    if not exists(out3_file):
        with open(out3_file, 'w', encoding='utf-8') as out3_title:
            out3_title.write(csv_title)
    out3 = open(out3_file, 'a', encoding='utf-8')
    out2.write(csv_title)
    for pair in pairs:
        line = style.format(
            locus, pair.score, rows, average(list(pair.length.values())),
            np.std(list(pair.length.values())), min(pair.length.values()),
            max(pair.length.values()),
            pair.coverage, pair.resolution, pair.tree_value,
            pair.avg_terminal_len, pair.entropy,
            pair.left.seq, pair.left.tm, pair.left.avg_bitscore,
            pair.left.avg_mismatch,
            pair.right.seq, pair.right.tm, pair.right.avg_bitscore,
            pair.right.avg_mismatch,
            pair.delta_tm, pair.left.start, pair.right.end, pair.start,
            pair.end)
        out2.write(line)
        out3.write(line)
        SeqIO.write(pair.left, out1, 'fastq')
        SeqIO.write(pair.right, out1, 'fastq')
    out1.close()
    out2.close()
    out3.close()
    log.info('Primers info were written into {}.csv.'.format(arg.out_file))
    return True


def analyze_wrapper(files, arg):
    """
    Wrapper for the primer design.
    """
    log.info('Analyze alignments.')
    result = []
    for aln in files:
        log.info('Analyze {}.'.format(aln))
        arg.out_file = splitext(clean_path(aln, arg))[0]
        result.append(analyze(aln, arg))
        log.info('')
    log.info('Analysis finished.')
    return any(result)


def quit(msg):
    """
    Quit for critical situation.
    """
    log.critical(msg)
    log.info('Quit.')
    raise SystemExit(-1)


def main():
    """
    main function
    """
    log.info('Welcome to BarcodeFinder v{}.'.format(__version__))
    arg = parse_args()
    if arg is None:
        quit('Empty input! Please use "-h" options for help info.')
    try:
        mkdir(arg.out)
    except FileExistsError:
        quit('Output folder "{}" is existed. Please use "-out" option to set '
             'a new output folder.'.format(arg.out))
    # move log to new log file
    log_tmp.close()
    log.removeHandler(log_tmp)
    log_file = join_path(arg.out, 'Log.txt')
    with open(TEMP_LOG, 'r') as a, open(log_file, 'w') as b:
        b.write(a.read())
    remove(TEMP_LOG)
    log_file_handler = logging.FileHandler(log_file)
    log_file_handler.setLevel(logging.INFO)
    log_file_handler.setFormatter(LOG_FMT)
    log.addHandler(log_file_handler)

    # prepare
    wrote_by_gene = []
    wrote_by_name = []
    mkdir(arg.by_gene_folder)
    mkdir(arg.by_name_folder)
    # collect and preprocess
    query = get_query_string(arg)
    if query is not None:
        log.info('Download data from Genbank.')
        gbfile = download(arg, query)
        if gbfile is not None:
            wrote_by_gene, wrote_by_name = divide(gbfile, arg)
        else:
            log.critical('Query is empty.')
    if arg.gb is not None:
        for i in list(glob(arg.gb)):
            by_gene, by_name = divide(i, arg)
            wrote_by_gene.extend(by_gene)
            wrote_by_name.extend(by_name)
    if arg.fasta is not None:
        user_data = list(glob(arg.fasta))
        wrote_by_name.extend(user_data)
    if not any([wrote_by_gene, wrote_by_name, arg.aln]):
        log.critical('Data is empty, please check your input!')
        log.info('Quit.')
        raise SystemExit(-1)
    if arg.uniq == 'no':
        log.info('Skip removing redundant sequences.')
    else:
        log.info('Remove redundant sequences by "{}".'.format(arg.uniq))
        wrote_by_gene = uniq(wrote_by_gene, arg)
        wrote_by_name = uniq(wrote_by_name, arg)
    if arg.stop == 1:
        log.info('Exit.')
        return
    # check dependent
    log.info('Checking dependent software.')
    original_path = check_tools()
    if original_path is None:
        quit('Cannot find and install depedent software.')
    # evaluate
    # only consider arg.no_divide and arg.fasta
    if arg.no_divide or arg.fasta:
        aligned = align(wrote_by_name, arg)
    else:
        aligned = align(wrote_by_gene, arg)
    # assume that alignments user provided is clean and do not nead uniq
    if arg.aln is not None:
        user_aln = list(glob(arg.aln))
        aligned.extend(user_aln)
    if len(aligned) == 0:
        log.critical('Cannot find valid alignment.')
        log.info('Quit')
        raise SystemExit(-1)
    result = analyze_wrapper(aligned, arg)
    log.info('Finished. You can find output in {}.'.format(arg.out))
    if result:
        log.info('Summary info were written into {} and {}.'.format(join_path(
            arg.out, 'Loci.csv'), join_path(arg.out, 'Primers.csv')))
    else:
        log.critical('None of input is valid.')
        log.info('Quit.')
        raise SystemExit(-1)
    json_file = join_path(arg.out, 'Options.json')
    with open(json_file, 'w', encoding='utf-8') as out:
        json.dump(vars(arg), out, indent=4, sort_keys=True)
    log.info('Options were dumped into {}.'.format(json_file))
    log.info('Reset PATH to original value.')
    environ['PATH'] = original_path
    log.info('Make sure temporary files were cleaned.')
    try:
        remove(arg.db_file)
    except FileNotFoundError:
        pass
    log.info('Bye.')
    return


__version__ = '0.9.35'
if __name__ == '__main__':
    main()
