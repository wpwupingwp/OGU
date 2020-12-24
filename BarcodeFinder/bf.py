#!/usr/bin/python3

import argparse
import logging
import json
import re
from collections import defaultdict
from glob import glob
from io import StringIO
from os import (cpu_count, devnull, environ, mkdir, pathsep, remove, rename,
                sep)
from os.path import abspath, basename, exists, splitext
from os.path import join as join_path
from random import choice
from subprocess import run
from time import sleep

import numpy as np
from primer3 import calcTm, calcHairpinTm, calcHomodimerTm, calcHeterodimerTm
from Bio import Entrez, Phylo, SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline as Blast
from Bio.Data.IUPACData import ambiguous_dna_values as ambiguous_data
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from BarcodeFinder import utils
from BarcodeFinder import gb2fasta
from BarcodeFinder import primer

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
DATEFMT = '%H:%M:%S'
logging.basicConfig(format=FMT, datefmt=DATEFMT, level=logging.INFO)
log = logging.getLogger('barcodefinder')
try:
    import coloredlogs
    coloredlogs.install(level=logging.INFO, fmt=FMT, datefmt=DATEFMT)
except ImportError:
    pass


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
            self.avg_mid_loc = int(utils.safe_average(list(
                self.mid_loc.values())))
        self.id = ('AvgMidLocation({:.0f})-Tm({:.2f})-Coverage({:.2%})-'
                   'AvgBitScore({:.2f})-Start({})-End({})'.format(
                    self.avg_mid_loc, self.tm, self.coverage,
                    self.avg_bitscore, self.start, self.end))





def parse_args():
    """
    Parse args and store some global/temporary values.
    """
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=main.__doc__)
    arg.add_argument('action',
                     choices=('all', 'gb2fasta', 'evaluate', 'primer'),
                     help=('all: do all analysis\t'
                           'gb2fasta: collect and organize\t'
                           'evaluate: align and evaluate\t'
                           'primer: design universal primer'))
    general = arg.add_argument_group('General')
    general.add_argument('-aln', help='aligned fasta files to analyze')
    general.add_argument('-fasta', help='unaligned fasta format data to add')
    general.add_argument('-gb', help='genbank files')
    general.add_argument('-out', help='output directory')
    evaluate = arg.add_argument_group('Evaluate')
    evaluate.add_argument('-fast', action='store_true', default=False,
                          help='faster evaluate variance by omit tree_value'
                          'and terminal branch length')
    evaluate.add_argument('-step', default=50, type=int,
                          help='step length for sliding-window scan')
    primer = arg.add_argument_group('Primer')
    primer.add_argument('-a', dest='ambiguous_base_n', default=4, type=int,
                        help='number of ambiguous bases')
    primer.add_argument('-c', dest='coverage', default=0.6, type=float,
                        help='minium coverage of base and primer')
    primer.add_argument('-m', dest='mismatch', default=4, type=int,
                        help='maximum mismatch bases in primer')
    primer.add_argument('-pmin', dest='min_primer', default=18, type=int,
                        help='minimum primer length')
    primer.add_argument('-pmax', dest='max_primer', default=28, type=int,
                        help='maximum primer length')
    primer.add_argument('-r', dest='resolution', type=float, default=0.5,
                        help='minium resolution')
    primer.add_argument('-t', dest='top_n', type=int, default=1,
                        help='keep n primers for each high varient region')
    primer.add_argument('-tmin', dest='min_product', default=350, type=int,
                        help='minimum product length(include primer)')
    primer.add_argument('-tmax', dest='max_product', default=600, type=int,
                        help='maximum product length(include primer)')
    parsed = arg.parse_args()
    if parsed.fast:
        log.info('The "-fast" mode was opened. '
                 'Skip sliding-window scan with tree.')
    # temporary filename, omit one parameters in many functions
    parsed.db_file = join_path(parsed.out, 'interleaved.fasta')
    parsed.out_file = ''
    # load option.json may cause chaos, remove
    return parsed


def clean_path(old, arg):
    """
    Join path if the file is not under by-gene or by-uniq to make working
    folder clean.
    """
    # to be continued
    try:
        split = old.split(sep)
    except:
        log.critical('Halt.')
        raise
    if 'by-gene' not in split and 'by-name' not in split:
        return join_path(arg.by_name_folder, basename(old))
    else:
        return old











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
        values = alignment.get_resolution(
            alignment, i, i + max_plus, arg.fast)
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
     a_branch_len) = alignment.get_resolution(alignment, 0, columns)
    log.info('\tGap ratio:\t{}'.format(a_gap_ratio))
    log.info('\tObserved resolution:\t{}'.format(a_observed_res))
    log.info('\tNormalized Shannon Index:\t{}'.format(a_entropy))
    log.info('\tPi:\t{}'.format(a_pi))
    log.info('\tTree resolution:\t{}'.format(a_tree_res))
    log.info('\tAverage terminal branch length:\t{}'.format(a_branch_len))
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
    good_region = primer.get_good_region(
        index, observed_res_list, arg)
    log.info('Finding candidate primers.')
    consensus = primer.find_continuous(
        consensus, good_region, arg.min_primer)
    primer_candidate, consensus = primer.find_primer(consensus, arg)
    if len(primer_candidate) == 0:
        log.warning('Cannot find primer candidates. '
                    'Please consider to loose options.')
        return True
    log.info('Found {} candidate primers.'.format(len(primer_candidate)))
    log.info('Validate with BLAST. May be slow.')
    primer_verified = primer.validate(primer_candidate, db_file, rows, arg)
    if len(primer_verified) == 0:
        log.warning('All candidates failed on validation. '
                    'Please consider to loose options.')
        return True
    log.info('Picking primer pairs.')
    pairs = primer.pick_pair(primer_verified, alignment, arg)
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
            locus, pair.score, rows, utils.safe_average(
                list(pair.length.values())),
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
    log.info('Welcome to BarcodeFinder.')
    arg = parse_args()
    init_ok, arg = init_arg(arg)
    log_file = arg.out / 'Log.txt'
    log_file_handler = logging.FileHandler(log_file, mode='a')
    log_file_handler.setLevel(logging.INFO)
    log.addHandler(log_file_handler)

    # prepare
    wrote_by_gene = []
    wrote_by_name = []
    mkdir(arg.by_gene_folder)
    mkdir(arg.by_name_folder)
    # collect and preprocess
    if arg.gb is not None:
        for i in list(glob(arg.gb)):
            by_gene, by_name = gb2fasta.divide(i, arg)
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
    if arg.stop == 1:
        log.info('Exit.')
        return
    # check dependent
    log.info('Checking dependent software.')
    original_path = utils.check_tools()
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
    log.info('Clean temporary files.')
    try:
        remove(arg.db_file)
    except FileNotFoundError:
        pass
    log.info('Bye.')
    return


if __name__ == '__main__':
    main()
