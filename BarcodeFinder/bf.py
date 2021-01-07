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
from BarcodeFinder import evaluate
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


def init_arg(arg):
    if arg.out is not None:
        arg.out = utils.init_out(arg)
    if not any([arg.gb, arg.fasta, arg.aln]):
        pass
    return arg


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

    if arg.action == 'all':
        pass
    elif arg.action == 'gb2fasta':
        gb2fasta.gb2fasta_main()
        return
    elif arg.action == 'evaluate':
        evaluate.evaluate_main()
        return
    elif arg.action == 'primer':
        primer.primer_main()
        return

    # collect and preprocess
    if not any([wrote_by_gene, wrote_by_name, arg.aln]):
        log.critical('Data is empty, please check your input!')
        log.info('Quit.')
        raise SystemExit(-1)
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
