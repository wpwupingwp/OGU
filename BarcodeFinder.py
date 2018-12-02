#!/usr/bin/python3

import argparse
import json
import re
from collections import defaultdict
from datetime import datetime
from glob import glob
from itertools import product as cartesian_product
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
    print('Cannot find DISPLAY. use Agg backend instead.')
    mpl_use('Agg')
from matplotlib import pyplot as plt
from matplotlib import rcParams

rcParams['axes.labelsize'] = 16
rcParams['axes.linewidth'] = 1.5
rcParams['axes.titlesize'] = 25
rcParams['font.size'] = 16
rcParams['lines.linewidth'] = 1.5


class PrimerWithInfo(SeqRecord):
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
        if isinstance(i, int):
            i = slice(i, i + 1)
        if isinstance(i, slice):
            answer = PrimerWithInfo(seq=str(self.seq[i]),
                                    quality=self.quality[i])
            answer.annotations = dict(self.annotations.items())
            return answer
        else:
            raise IndexError

    def is_good_primer(self):
        poly = re.compile(r'([ATCG])\1\1\1\1')
        tandem = re.compile(r'([ATCG]{2})\1\1\1\1')
        # ref1. http://www.premierbiosoft.com/tech_notes/PCR_Primer_Design.html
        if re.search(poly, str(self.seq)) is not None:
            self.detail = 'Poly(NNNNN) structure found'
            return False
        if re.search(tandem, str(self.seq)) is not None:
            self.detail = 'Tandom(NN*5) exist'
            return False
        self.hairpin_tm = calc_ambiguous_seq(calcHairpinTm, self.seq)
        self.homodimer_tm = calc_ambiguous_seq(calcHomodimerTm, self.seq)
        self.tm = self.annotations['tm'] = calc_ambiguous_seq(calcTm, self.seq)
        # primer3.calcHairpin or calcHomodimer usually return structure found
        # with low Tm. Here we compare structure_tm with sequence tm
        if self.hairpin_tm >= self.tm:
            self.detail = 'Hairpin found'
            return False
        if self.homodimer_tm >= self.tm:
            self.detail = 'Homodimer found'
            return False
        return True

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
                 'have_heterodimer', 'heterodimer_tm', 'pi', 'score', 'length']

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
        self.score = 0.0
        self.get_score()

    def __str__(self):
        return (
            'Pair(score={:.2f}, product={:.0f}, start={}, end={}, left={}, '
            'right={}, observerd_resolution={:.2%}, coverage={:.2%},'
            'delta_tm={:.2f}, have_heterodimer={})'.format(
                self.score, average(list(self.length.values())), self.start,
                self.end, self.left.seq, self.right.seq, self.resolution,
                self.coverage, self.delta_tm, self.have_heterodimer))

    def get_score(self):
        self.score = (average(list(self.length.values())) * 0.5
                      + self.coverage * 200
                      + len(self.left) * 10
                      + len(self.right) * 10
                      + self.resolution * 100
                      + self.tree_value * 100 + self.entropy * 5
                      - int(self.have_heterodimer) * 10
                      - self.delta_tm * 5 - self.left.avg_mismatch * 10
                      - self.right.avg_mismatch * 10)

    def add_info(self, alignment):
        """
        put slow steps here to save time
        """
        if not self.right.is_reverse_complement:
            self.right = self.right.reverse_complement()
        # include end base, use alignment loc for slice
        (self.resolution, self.entropy, self.pi,
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
    __slots = ('query_id', 'hit_id', 'query_seq', 'ident_num', 'mismatch_num',
               'bitscore_raw', 'query_start', 'query_end', 'hit_start',
               'hit_end')

    def __init__(self, line):
        record = line.strip().split('\t')
        self.query_id, self.hit_id, self.query_seq = record[0:3]
        (self.ident_num, self.mismatch_num, self.bitscore_raw,
         self.query_start, self.query_end, self.hit_start,
         self.hit_end) = [int(i) for i in record[3:]]


def parse_args():
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=main.__doc__)
    general = arg.add_argument_group('General')
    general.add_argument('-aln', help='aligned fasta files to analyze')
    general.add_argument('-fasta', help='unaligned fasta format data to add')
    general.add_argument('-gb', help='genbank files')
    general.add_argument('-stop', type=int, choices=(1, 2, 3), default=3,
                         help=('Stop after which step:'
                               '\t1. Download and pre-process;'
                               '\t2. Analyze variance;'
                               '\t3. Primer design.'))
    general.add_argument('-out', help='output directory')
    genbank = arg.add_argument_group('Genbank')
    genbank.add_argument('-email', help='email address for querying Genbank')
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
    genbank.add_argument('-organelle',
                         choices=('mitochondrion', 'plastid', 'chloroplast'),
                         help='organelle type')
    genbank.add_argument('-query', help='query text')
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
    evaluate.add_argument('-f', dest='fast', action='store_true',
                          default=False,
                          help='faster evaluate variance by omit tree_value'
                          'and terminal branch length')
    evaluate.add_argument('-s', dest='step', type=int, default=50,
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
    if parsed.organelle is not None:
        # 10k to 1m seems enough
        parsed.min_len = 10000
        parsed.max_len = 1000000
    if not any([parsed.query, parsed.taxon, parsed.gene, parsed.fasta,
                parsed.aln, parsed.gb]):
        arg.print_help()
        raise ValueError('Empty input!')
    if parsed.out is None:
        raw_time = datetime.now().isoformat()
        parsed.out = raw_time.replace(':', '-').split('.')[0]
    parsed.by_gene_folder = join_path(parsed.out, 'by-gene')
    parsed.by_name_folder = join_path(parsed.out, 'by-name')
    # temporary filename, omit one parameters in many functions
    parsed.db_file = join_path(parsed.out, 'interleaved.fasta')
    parsed.no_gap_file = join_path(parsed.out, 'no_gap.fasta')
    parsed.out_file = ''
    # load option.json may cause chaos, remove
    return parsed


def tprint(string):
    now = datetime.now()
    s = '{:0>2d}:{:0>2d}:{:>02d}   {}'.format(now.hour, now.minute, now.second,
                                              string)
    print(s, flush=True)
    log_handle.write(s + '\n')
    log_handle.flush()


def average(x):
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
    Join path if the file is not under by-gene or by-uniq.
    To make working folder clean.
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
    # Seems primer3 only accept seqs shorter than 60 bp.
    # Plus, too long seq will cost too much memory
    len_limit = 60

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

    if len(seq) > len_limit:
        return 0
    seq_str = _expand(seq)
    if seq2 is None:
        values = [func(i) for i in seq_str]
    else:
        if len(seq2) > len_limit:
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
            tprint('Cannot find {}. Try to install.'.format(tools))
            install_path = deploy(tools)
            if install_path is None:
                tprint('Failed to install {}. Please try to manually install'
                       'it (See README.md).'.format(tools))
                return None
            installed.append(install_path)
    # do not edit original PATH
    to_add = pathsep.join(installed)
    original = str(environ['PATH'])
    environ['PATH'] = pathsep.join([original, to_add])
    with open('PATH.txt', 'w', encoding='utf-8') as path_out:
        path_out.write(to_add + '\n')
    f.close()
    return original


def download_software(url):
    """
    Download, return False if failed.
    http_proxy may affect this function.
    """
    filename = url.split('/')[-1]
    try:
        down = urlopen(url)
    except HTTPError:
        tprint('Cannot download {}.'.format(filename))
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
    return False if failed
    """
    tprint('Try to install {}. Please consider to install it following '
           'official instruction to get a CLEAN system.'.format(software))
    sys = system()
    # url dict
    blast_url = ('ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/'
                 'ncbi-blast-2.7.1+')
    iqtree_url = ('https://github.com/Cibiv/IQ-TREE/releases/download/v1.6.8/'
                  'iqtree-1.6.8')
    mafft_url = 'https://mafft.cbrc.jp/alignment/software/mafft'
    # windows blast path not sure
    urls = {'Linux':
            {'BLAST': {'url': blast_url+'-x64-linux.tar.gz',
                       'path': abspath('ncbi-blast-2.7.1+'+sep+'bin')},
             'IQTREE': {'url': iqtree_url+'-Linux.tar.gz',
                        'path': abspath('iqtree-1.6.8-Linux'+sep+'bin')},
             'MAFFT': {'url': mafft_url+'-7.407-linux.tgz',
                       'path': abspath('mafft-linux64')}},
            'macOSX':
            {'BLAST': {'url': blast_url+'.dmg',
                       'path': abspath('ncbi-blast-2.7.1+'+sep+'bin')},
             'IQTREE': {'url': iqtree_url+'-MacOSX.zip',
                        'path': abspath('iqtree-1.6.8-MacOSX'+sep+'bin')},
             'MAFFT': {'url': mafft_url+'-7.407-mac.zip',
                       'path': abspath('mafft-mac')}},
            'Windows':
            {'BLAST': {'url': blast_url+'-win64.exe',
                       'path': abspath('.')},
             'IQTREE': {'url': iqtree_url+'-Windows.zip',
                        'path': abspath('iqtree-1.6.8-Windows'+sep+'bin')},
             'MAFFT': {'url': mafft_url+'-7.409-win64-signed.zip',
                       'path': abspath('mafft-win')}}}
    url = urls[sys][software]['url']
    # down
    if sys == 'Windows':
        if not download_software(url):
            return None
        if software == 'BLAST':
            run('ncbi-blast-2.7.1+-win64.exe', shell=True)
    elif sys == 'Linux':
        ok = False
        for pack_mgr in ('apt', 'dnf', 'yum', 'pkg'):
            r = run('sudo {} install ncbi-blast+ iqtree mafft'.format(
                pack_mgr), shell=True)
            if r.returncode == 0:
                ok = True
                break
        if not ok:
            tprint('Cannot install {} to system, try to'
                   'download.'.format(software))
            download_software(url)
    elif sys == 'macOSX':
        brew_ok = False
        ok = False
        with open(devnull, 'w', encoding='utf-8') as f:
            r = run('brew --help', shell=True, stdout=f, stderr=f)
        if r.returncode == 0:
            brew_ok = True
        else:
            r2 = run('/usr/bin/ruby -e "$(curl -fsSL https://raw.'
                     'githubusercontent.com/Homebrew/install/master/install)"',
                     shell=True)
            if r2.returncode == 0:
                brew_ok = True
            else:
                tprint('Cannot install brew.')
        if brew_ok:
            r3 = run('brew install blast mafft brewsci/science/iqtree',
                     shell=True)
            if r3.returncode == 0:
                ok = True
        if not ok:
                download_software(url)
    # windows can omit .bat, linux cannot
    if software == 'MAFFT' and sys != 'Windows':
        rename(join_path(urls[sys]['mafft']['path'], 'mafft.bat'),
               join_path(urls[sys]['mafft']['path'], 'mafft'))
    return abspath(urls[sys][software]['path'])


def get_query_string(arg):
    condition = []
    if arg.group is not None:
        condition.append('{}[filter]'.format(arg.group))
    if arg.query is not None:
        condition.append('"{}"'.format(arg.query))
    if arg.gene is not None:
        condition.append('"{}"[gene]'.format(arg.gene))
    if arg.molecular is not None:
        d = {'DNA': 'biomol_genomic[PROP]',
             'RNA': 'biomol_mrna[PROP]'}
        condition.append(d[arg.molecular])
    if arg.taxon is not None:
        condition.append('"{}"[ORGANISM]'.format(arg.taxon))
    if arg.organelle is not None:
        condition.append('{}[filter]'.format(arg.organelle))
    if len(condition) != 0:
        condition.append('("{}"[SLEN] : "{}"[SLEN])'.format(
            arg.min_len, arg.max_len))
        return ' AND '.join(condition)
    else:
        return None


def download(arg, query):
    tprint('Query:\t{}.'.format(query))
    if arg.email is None:
        Entrez.email = 'guest@example.com'
        tprint('You did not provide email address, use'
               ' {} instead.'.format(Entrez.email))
    else:
        Entrez.email = arg.email
    query_handle = Entrez.read(Entrez.esearch(db='nuccore', term=query,
                                              usehistory='y'))
    count = int(query_handle['Count'])
    tprint('{} records found.'.format(count))
    tprint('Downloading... Ctrl+C to quit.')
    json_file = join_path(arg.out, 'Query.json')
    with open(json_file, 'w', encoding='utf-8') as _:
        json.dump(query_handle, _, indent=4, sort_keys=True)
    if arg.query is None:
        name = safe(arg.taxon)
    else:
        name = safe(arg.query)
    file_name = join_path(arg.out, name + '.gb')
    output = open(file_name, 'w', encoding='utf-8')
    ret_start = 0
    if count >= 1000:
        ret_max = 1000
    elif count >= 100:
        ret_max = 100
    else:
        ret_max = 10
    while ret_start < count:
        tprint('{:d}--{:d}'.format(ret_start, ret_start + ret_max))
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
            tprint('Retrying...')
            continue
        ret_start += ret_max
    tprint('Download finished.')
    return file_name


def gene_rename(old_name):
    """For chloroplast gene.
    Input->str
    Output->List[new_name:str, name_type:str]
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


def get_taxon(order_family, order_exceptions):
    # order|family|organims(genus|species)
    order = ''
    family = ''
    for item in order_family:
        if item.endswith('ales') or item in order_exceptions:
            order = item
        elif item.endswith('aceae') or item.endswith('idae'):
            family = item
    return order, family


def write_seq(name, sequence_id, feature, whole_seq, path, arg):
    """
    Write fasta file.
    """

    def careful_extract(whole_seq):
        try:
            sequence = feature.extract(whole_seq)
        except ValueError:
            sequence = ''
            tprint('Cannot extract sequence of {} from {}.'.format(
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
    Return: [name, feature.type]
    """
    name = None
    misc_feature = None
    if feature.type == 'gene':
        if 'gene' in feature.qualifiers:
            gene = feature.qualifiers['gene'][0].replace(' ', '_')
            if arg.rename:
                gene = gene_rename(gene)[0]
            name = safe(gene)
        elif 'product' in feature.qualifiers:
            product = feature.qualifiers['product'][0].replace(
                ' ', '_')
            name = safe(product)
    elif feature.type == 'misc_feature':
        if 'product' in feature.qualifiers:
            misc_feature = feature.qualifiers['product'][0].replace(
                ' ', '_')
        elif 'note' in feature.qualifiers:
            misc_feature = feature.qualifiers['note'][0].replace(
                ' ', '_')
        if (misc_feature is not None) and ('intergenic_spacer' in misc_feature
                                           or 'IGS' in misc_feature):
            # 'IGS' in misc_feature) and len(misc_feature) < 100):
            name = safe(misc_feature)
            name = name.replace('intergenic_spacer_region',
                                'IGS')
    elif feature.type == 'misc_RNA':
        if 'product' in feature.qualifiers:
            misc_feature = feature.qualifiers['product'][0].replace(
                ' ', '_')
        elif 'note' in feature.qualifiers:
            misc_feature = feature.qualifiers['note'][0].replace(
                ' ', '_')
        name = safe(misc_feature)
        # handle ITS
        if 'internal_transcribed_spacer' in name:
            name = 'ITS'
        # name = name.replace('internal_transcribed_spacer', 'ITS')
        # if 'ITS_1' in name:
        #     if 'ITS_2' in name:
        #         name = 'ITS'
        #     else:
        #         name = 'ITS_1'
        # elif 'ITS_2' in name:
        #     name = 'ITS_2'
    else:
        pass
    return name, feature.type


def get_spacer(genes):
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
    # put raw fasta into root of output folder, so not to use clean_path
    raw_fasta = join_path(arg.out, splitext(basename(gbfile))[0] + '.fasta')
    handle_raw = open(raw_fasta, 'w', encoding='utf-8')
    wrote_by_gene = set()
    wrote_by_name = set()
    # order exceptions that not end with "ales"
    # From NCBI Taxonomy Database, 2018.11.27 update
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
    order_exceptions = order_exceptions.split(',')
    order_exceptions = [i.strip() for i in order_exceptions]
    # divide gb
    for record in SeqIO.parse(gbfile, 'gb'):
        # only accept gene, product, and spacer in misc_features.note
        order_family = record.annotations['taxonomy']
        order, family = get_taxon(order_family, order_exceptions)
        organism = record.annotations['organism'].replace(' ', '_')
        genus, *species = organism.split('_')
        # species name may contain other characters
        taxon = '{}|{}|{}|{}'.format(order, family, genus, '_'.join(species))
        accession = record.annotations['accessions'][0]
        try:
            specimen = record.features[0].qualifiers['specimen_voucher'
                                                     ][0].replace(' ', '_')
        except (IndexError, KeyError):
            specimen = ''
        whole_seq = record.seq
        feature_name = []
        genes = []

        for feature in record.features:
            name, feature_type = get_feature_name(feature, arg)
            # skip unsupport feature
            if name is None:
                continue
            if len(name) > arg.max_name_len:
                tprint('Too long name: {}.'.format(name))
                name = name[:arg.max_name_len] + '...'
            # skip abnormal annotation
            if len(feature) > arg.max_seq_len:
                tprint('Skip abnormal annotaion of {}(Accession {}).'.format(
                    name, accession))
                continue
            if feature_type == 'gene':
                genes.append([name, feature])
            feature_name.append(name)
            sequence_id = '>' + '|'.join([name, taxon, accession, specimen])
            wrote = write_seq(name, sequence_id, feature, whole_seq,
                              arg.by_gene_folder, arg)
            wrote_by_gene.add(wrote)

        # extract spacer
        spacers = get_spacer(genes)
        for spacer in spacers:
            if len(spacer) > arg.max_seq_len:
                tprint('Spacer {} too long (Accession {}).'.format(
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
        wrote_by_name.remove(unknown)
    tprint('Divide done.')
    return list(wrote_by_gene), list(wrote_by_name)


def uniq(files, arg):
    uniq_files = []
    for fasta in files:
        info = defaultdict(lambda: list())
        keep = dict()
        count = 0
        for record in SeqIO.parse(fasta, 'fasta'):
            # gene|order|family|genus|species|specimen
            if '|' in record.id:
                name = ' '.join(record.id.split('|')[3:5])
            else:
                name = record.id
            length = len(record)
            # skip empty file
            if length != 0:
                info[name].append([count, length])
            count += 1
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
        elif arg.uniq == 'no':
            # keep all
            keep = {range(count + 1)}
        new = clean_path(fasta, arg) + '.uniq'
        with open(new, 'w', encoding='utf-8') as out:
            for index, record in enumerate(SeqIO.parse(fasta, 'fasta')):
                if index in keep:
                    SeqIO.write(record, out, 'fasta')
        uniq_files.append(new)
    return uniq_files


def align(files, arg):
    result = []
    # get available CPU cores
    cores = cpu_count() - 1
    for fasta in files:
        tprint('Aligning {}.'.format(fasta))
        out = clean_path(fasta, arg) + '.aln'
        with open(devnull, 'w', encoding='utf-8') as f:
            _ = ('mafft --thread {} --reorder --quiet --adjustdirection '
                 '{} > {}'.format(cores, fasta, out))
            m = run(_, shell=True, stdout=f, stderr=f)
        if m.returncode == 0:
            result.append(out)
        else:
            tprint('Skip alignment of {}.'.format(fasta))
    tprint('Alignment done.')
    for i in glob('_order*'):
        remove(i)
    return result


def prepare(aln_fasta, arg):
    """
    Given fasta format alignment filename, return a numpy array for sequence:
    Generate fasta file without gap for makeblastdb, return file name.
    """
    data = []
    record = ['id', 'sequence']
    with open(aln_fasta, 'r', encoding='utf-8') as raw, open(
            arg.no_gap_file, 'w', encoding='utf-8') as no_gap:
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
        tprint('{} does not have uniform width!'.format(aln_fasta))
        return None, None, None

    # Convert List to numpy array.
    # order 'F' is a bit faster than 'C'
    # new = np.hstack((name, seq)) -> is slower
    name = np.array([[i[0]] for i in data], dtype=np.bytes_)
    sequence = np.array([list(i[1]) for i in data], dtype=np.bytes_, order='F')

    if name is None:
        tprint('Bad fasta file {}.'.format(aln_fasta))
        name = None
    # tree require more than 4 sequences
    if len(sequence) < 4:
        tprint('Too few sequence in {} (less than 4)!'.format(aln_fasta))
        name = None
    interleaved = 'interleaved.fasta'
    # for clean
    # try to avoid makeblastdb error
    SeqIO.convert(arg.no_gap_file, 'fasta', interleaved, 'fasta')
    return name, sequence, interleaved


def count_base(alignment, rows, columns):
    """
    Given alignment numpy array, count cumulative frequency of base in each
    column (consider ambiguous base and "N", "-" and "?", otherwise omit).
    Return List[List[float, float, float, float, float, float, float]] for
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
    return resolution, entropy, Pi and tree value.
    """
    subalign = alignment[:, start:end]
    rows, columns = subalign.shape
    # index error
    if columns == 0:
        return 0, 0, 0, 0
    item, count = np.unique(subalign, return_counts=True, axis=0)
    resolution = len(count) / rows
    tree_value = 0
    avg_terminal_branch_len = 0
    # entropy
    entropy = 0
    for j in count:
        p_j = j / rows
        log2_p_j = np.log2(p_j)
        entropy += log2_p_j * p_j
    entropy *= -1
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
        with open(aln_file, 'wb') as aln:
            for index, row in enumerate(alignment[:, start:end]):
                aln.write(b'>' + str(index).encode('utf-8') + b'\n' + b''.join(
                    row) + b'\n')
        with open(devnull, 'w', encoding='utf-8') as f:
            iqtree = run('iqtree -s {} -m JC -fast -czb'.format(aln_file),
                         stdout=f, stderr=f, shell=True)
        # just return 0 if there is error
        if iqtree.returncode != 0:
            tprint('Too much gap in the region {}-{} bp.'.format(
                start, end))
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
    return resolution, entropy, pi, tree_value, avg_terminal_branch_len


def generate_consensus(base_cumulative_frequency, coverage_percent,
                       rows, output):
    """
    Given base count info, return List[index, base, quality]
    and List[List[str, str, str, PrimerInfo]] for writing consensus.
    return PrimerWithInfo
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

    limit = coverage / 4
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
    # return loose region, final product may violate product length
    # restriction
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
    Given PrimerWithInfo, good_region: List[bool], min_len
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


def find_primer(consensus, min_len, max_len):
    """
    Find suitable primer in given consensus with features labeled as candidate
    primer, return List[PrimerWithInfo], consensus
    """
    primers = []
    # skip good_region
    continuous = consensus.features
    for feature in continuous:
        fragment = feature.extract(consensus)
        len_fragment = len(fragment)
        for begin in range(len_fragment - max_len):
            for p_len in range(min_len, max_len + 1):
                start = feature.location.start + begin
                primer = consensus[start:start + p_len]
                if primer.is_good_primer():
                    consensus.features.append(SeqFeature(
                        FeatureLocation(start, start + p_len),
                        type='primer', strand=1))
                    primer.start = start
                    primer.update_id()
                    primers.append(primer)
    return primers, consensus


def count_and_draw(alignment, arg):
    """
    Given alignment(numpy array), return unique sequence count List[float].
    Calculate Shannon Index based on
    www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/shannon.htm
    return List[float]
    All calculation excludes primer sequence.
    """
    output = join_path(arg.out, basename(arg.out_file).split('.')[0])
    rows, columns = alignment.shape
    min_primer = arg.min_primer
    max_product = arg.max_product
    step = arg.step
    # r_list, h_list, pi_list, t_list : count, normalized entropy, Pi and
    #  tree value
    r_list = []
    h_list = []
    pi_list = []
    t_list = []
    l_list = []
    max_h = np.log2(rows)
    index = []
    max_plus = max_product - min_primer * 2
    max_range = columns - max_product
    for i in range(0, max_range, step):
        # skip gap
        # if consensus.sequence[i] in ('-', 'N'):
        #     continue
        # exclude primer sequence
        resolution, entropy, pi, tree_value, avg_t_branch_len = get_resolution(
            alignment, i, i + max_plus, arg.fast)
        r_list.append(resolution)
        h_list.append(entropy / max_h)
        pi_list.append(pi)
        t_list.append(tree_value)
        l_list.append(avg_t_branch_len)
        index.append(i)

    plt.style.use('seaborn-colorblind')
    # how to find optimized size?
    fig, ax1 = plt.subplots(figsize=(15 + len(index) // 5000, 10))
    plt.title('Variance of {} (window={} bp, step={} bp)\n'.format(
        basename(arg.out_file).split('.')[0], max_product, step))
    plt.xlabel('Base')
    ax1.yaxis.set_ticks(np.linspace(0, 1, num=11))
    ax2 = ax1.twinx()
    ax2.set_ylabel(r'$\pi$', rotation=-90, labelpad=20)
    # plt.xticks(np.linspace(0, max_range, 21))
    if not arg.fast:
        ax1.plot(index, t_list, label='Tree Resolution', alpha=0.8)
        ax2.plot(index, l_list, linestyle='--',
                 label='Average Terminal Branch Length', alpha=0.8)
        ax2.set_ylabel(r'$\pi$ and Average Branch Length', rotation=-90,
                       labelpad=20)
    ax1.set_ylabel('Resolution')
    ax1.plot(index, h_list, label='Shannon Equitability Index', alpha=0.8)
    ax1.plot(index, r_list, label='Observed Resolution', alpha=0.8)
    ax1.legend(loc='lower left')
    ax2.plot(index, pi_list, 'k-', label=r'$\pi$', alpha=0.8)
    # _ = round(np.log10(max(Pi))) + 1
    # ax2.yaxis.set_ticks(np.linspace(0, 10**_, num=11))
    ax2.legend(loc='upper right')
    # plt.yscale('log')
    plt.savefig(output + '.pdf')
    plt.savefig(output + '.png')
    plt.close()
    # plt.show()
    with open(output + '.variance.tsv', 'w', encoding='utf-8') as _:
        _.write('Index,Resolution,TreeValue,AvgTerminalBranchLen,Entropy,Pi\n')
        for i, r, t, l, h, pi in zip(index, r_list, t_list, l_list, h_list,
                                     pi_list):
            # iqtree blmin is 1e-6
            _.write('{},{:.6f},{:.6f},{:.6f},{:.6f},{:.6f}\n'.format(i, r, t,
                                                                     l, h, pi))
    return r_list, h_list, pi_list, t_list, l_list, index


def parse_blast_tab(filename):
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
    Do BLAST. Parse BLAST result. Return List[PrimerWithInfo]
    """
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
            tprint('Failed to run makeblastdb!')
            return []
    # blast
    tprint('Validate with BLAST.')
    blast_result_file = 'blast.result.tsv'
    fmt = 'qseqid sseqid qseq nident mismatch score qstart qend sstart send'
    cmd = Blast(num_threads=cpu_count() - 1,
                query=query_file,
                db=db_file,
                task='blastn-short',
                evalue=1e-2,
                max_hsps=1,
                outfmt='"7 {}"'.format(fmt),
                out=blast_result_file)
    # hide output
    cmd()
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
    return primer_verified


def pick_pair(primers, alignment, arg):
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
    tprint('{} pairs of redundant primers were removed.'.format(
        len(pairs) - len(less_pairs)))
    good_pairs = []
    for i in less_pairs:
        i.add_info(alignment)
        if i.resolution >= arg.resolution:
            good_pairs.append(i)
    good_pairs.sort(key=lambda x: x.score, reverse=True)
    tprint('Successfully found validated primers.')
    return good_pairs[:arg.top_n]


def analyze(fasta, arg):
    """
    Automatic design primer for DNA barcode.
    """
    # read from fasta, generate new fasta for makeblastdb
    name, alignment, db_file = prepare(fasta, arg)
    if name is None:
        tprint('Invalid fasta file {}.'.format(fasta))
        return False
    rows, columns = alignment.shape
    # generate consensus
    base_cumulative_frequency = count_base(alignment, rows, columns)
    tprint('Generate consensus.')
    consensus = generate_consensus(base_cumulative_frequency, arg.coverage,
                                   rows, arg.out_file + '.consensus.fastq')
    tprint('Evaluate whole alignment.')
    max_count, max_h, max_pi, max_t, max_l = get_resolution(alignment, 0,
                                                            columns)
    tprint('Average terminal branch length {}.'.format(max_l))
    n_gap = sum([i[5] for i in base_cumulative_frequency])
    gap_ratio = n_gap / rows / columns
    summary = join_path(arg.out, 'Loci.csv')
    if not exists(summary):
        with open(summary, 'w', encoding='utf-8') as s:
            s.write('Loci,Sequences,Length,GapRatio,ObservedResolution,'
                    'TreeResolution,ShannonIndex,AvgTerminalBranchLen,Pi\n')
            s.write('{},{},{},{:.2%},{:.6f},{:.6f},{:.6f},{:.6f},{:.6f}'
                    '\n'.format(basename(fasta), rows, columns, gap_ratio,
                                max_count, max_t, max_h, max_l, max_pi))
    else:
        with open(summary, 'a', encoding='utf-8') as s:
            s.write('{},{},{},{:.2%},{:.6f},{:.6f},{:.6f},{:.6f},{:.6f}'
                    '\n'.format(basename(fasta), rows, columns, gap_ratio,
                                max_count, max_t, max_h, max_l, max_pi))
    if max_count < arg.resolution:
        tprint('Too low resolution of {}!'.format(fasta))
        return False
    # count resolution
    tprint('Sliding window analyze.')
    (seq_count, H, Pi, T, L, index) = count_and_draw(alignment, arg)
    # exit if resolution lower than given threshold.
    if len(seq_count) == 0:
        tprint('Problematic Input of {}!'.format(fasta))
    # stop if do not want to design primer
    if arg.stop == 2:
        return True
    # find ncandidate
    tprint('Start finding primers of {}.'.format(fasta))
    good_region = get_good_region(index, seq_count, arg)
    consensus = find_continuous(consensus, good_region, arg.min_primer)
    tprint('Filtering candidate primer pairs.')
    primer_candidate, consensus = find_primer(consensus, arg.min_primer,
                                              arg.max_primer)
    if len(primer_candidate) == 0:
        tprint('Cannot find primer candidates in {}. Try to loose'
               'options!'.format(fasta))
        return True
    tprint('Found {} candidate primers.'.format(len(primer_candidate)))
    # validate
    primer_verified = validate(primer_candidate, db_file, rows, arg)
    if len(primer_verified) == 0:
        tprint('Cannot find primers in {}. Try to loose options!'.format(
            fasta))
        return True
    # pick pair
    pairs = pick_pair(primer_verified, alignment, arg)
    if len(pairs) == 0:
        tprint('Cannot find primers in {}. Try to loose options!'.format(
            fasta))
        return True
    # output
    locus = basename(arg.out_file).split('.')[0]
    csv_title = ('Locus,Score,Sequences,AvgProductLength,StdEV,'
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
    tprint('Primers info were written into {}.csv.'.format(arg.out_file))
    return True


def analyze_wrapper(files, arg):
    result = []
    for aln in files:
        tprint('Analyze {}.'.format(aln))
        arg.out_file = splitext(clean_path(aln, arg))[0]
        result.append(analyze(aln, arg))
    # dirty work
    try:
        remove(arg.no_gap_file)
        remove(arg.db_file)
    except FileNotFoundError:
        pass
    return any(result)


def main():
    # prepare
    arg = parse_args()
    mkdir(arg.out)
    wrote_by_gene = []
    wrote_by_name = []
    mkdir(arg.by_gene_folder)
    mkdir(arg.by_name_folder)
    global log_handle
    log_handle = open(join_path(arg.out, 'Log.txt'), 'w', encoding='utf-8')
    tprint('Welcome to BarcodeFinder!')
    original_path = check_tools()
    if original_path is None:
        tprint('Exit.')
        return
    # collect and preprocess
    query = get_query_string(arg)
    if query is not None:
        tprint('Download data from Genbank.')
        gbfile = download(arg, query)
        tprint('Divide data by annotation.')
        wrote_by_gene, wrote_by_name = divide(gbfile, arg)
    if arg.gb is not None:
        for i in list(glob(arg.gb)):
            tprint('Divide {}.'.format(i))
            by_gene, by_name = divide(i, arg)
            wrote_by_gene.extend(by_gene)
            wrote_by_name.extend(by_name)
    if arg.fasta is not None:
        user_data = list(glob(arg.fasta))
        wrote_by_name.extend(user_data)
    if not any([wrote_by_gene, wrote_by_name, arg.aln]):
        tprint('Data is empty, please check your input!')
        environ['PATH'] = original_path
        return
    if arg.uniq == 'no':
        tprint('Skip removing redundant sequences.')
    else:
        tprint('Remove redundant sequences by "{}".'.format(arg.uniq))
    wrote_by_gene = uniq(wrote_by_gene, arg)
    wrote_by_name = uniq(wrote_by_name, arg)
    if arg.stop == 1:
        environ['PATH'] = original_path
        return
    # evaluate
    tprint('Aligning sequences.')
    # only consider arg.no_divide and arg.fasta
    if arg.no_divide or arg.fasta:
        aligned = align(wrote_by_name, arg)
    else:
        aligned = align(wrote_by_gene, arg)
    # assume that alignments user provided is clean and do not nead uniq
    if arg.aln is not None:
        user_aln = list(glob(arg.aln))
        aligned.extend(user_aln)
    result = analyze_wrapper(aligned, arg)
    tprint('Finished. You can find output in {}.'.format(arg.out))
    if result:
        tprint('Summary info were written into {} and {}.'.format(join_path(
            arg.out, 'Loci.csv'), join_path(arg.out, 'Primers.csv')))
    json_file = join_path(arg.out, 'Options.json')
    with open(json_file, 'w', encoding='utf-8') as out:
        json.dump(vars(arg), out, indent=4, sort_keys=True)
    tprint('Options were dumped into {}.'.format(json_file))
    log_handle.close()
    # restore original PATH
    environ['PATH'] = original_path
    return


if __name__ == '__main__':
    main()
