#!/usr/bin/python3

import argparse
import logging
import sys
import numpy as np
from pathlib import Path
from glob import glob
from os import remove, devnull, cpu_count
from subprocess import run

from Bio import Phylo

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


def parse_args(arg_list=None):
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=evaluate_main.__doc__)
    arg.add_argument('-fasta', nargs='*', help='unaligned fasta files')
    arg.add_argument('-aln', nargs='*', help='aligned files')
    arg.add_argument('-out', help='output folder')
    options = arg.add_argument_group('Options')
    options.add_argument('-ig', '-ignore_gap', action='store_true',
                         help='ignore gaps in alignments')
    options.add_argument('-iab', '-ignore_ambiguous_base', action='store_true',
                         help='ignore ambiguous bases like "M" or "N"')
    sliding_window = arg.add_argument_group('Sliding-window')
    sliding_window.add_argument('-fast', action='store_true',
                                help='skip evaluate phylogenetic diversity')
    sliding_window.add_argument('-size', type=int, default=500,
                                help='window size')
    sliding_window.add_argument('-step', type=int, default=50,
                                help='step length')
    if arg_list is None:
        return arg.parse_args()
    else:
        return arg.parse_args(arg_list)


def init_arg(arg):
    if arg.fast:
        log.info('The "-fast" mode was opened. '
                 'Skip evaluating phylogenetic diversity')
    return arg


def align(files: list, arg) -> (list, list):
    """
    Align sequences with mafft.
    Args:
        files(list): fasta files
        arg: arg
    Returns:
        aligned(list): aligned fasta
        unaligned(list): unaligned files

    """
    log.info('Align sequences.')
    aligned = list()
    unaligned = list()
    # get available CPU cores
    cores = max(1, cpu_count() - 1)
    for fasta in files:
        log.info('Aligning {}.'.format(fasta))
        out = arg._align / fasta.with_suffix('.aln').name
        with open(devnull, 'w', encoding='utf-8') as f:
            # if computer is good enough, "--genafpair" is recommended
            # where is mafft?
            _ = (f'{mafft} --auto --thread {} --reorder --quiet '
                 f'--adjustdirection {fasta} > {out}')
            m = run(_, shell=True, stdout=f, stderr=f)
        if m.returncode == 0:
            aligned.append(out)
        else:
            unaligned.append(fasta)
    log.info('Alignment finished.')
    for i in Path().cwd().glob('_order*'):
        i.unlink()
    log.info(f'Total {len(files)} files')
    log.info(f'Aligned {len(aligned)}')
    log.info(f'Unaligned {len(unaligned)}')
    return aligned, unaligned


def remove_gap(alignment):
    new_alignment = None
    return new_alignment


def get_resolution(alignment, start, end, fast=False):
    """
    Given alignment (2d numpy array), location of fragment(start and end, int,
    start from zero, exclude end),
    return gap ratio, resolution, entropy, Pi, tree value and average terminal
    branch length.
    """
    subalign = alignment[:, start:end]
    rows, columns = subalign.shape
    old_max_recursion = sys.getrecursionlimit()
    sys.setrecursionlimit(max(rows+10, old_max_recursion))
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
            try:
                internals = tree.get_nonterminals()[1:]
                terminals = tree.get_terminals()
                sum_terminal_branch_len = sum([i.branch_length for i in terminals])
            except RecursionError:
                log.info('Bad phylogentic tree.')
                internals = 0
                terminals = 1
                sum_terminal_branch_len = 0
            # miss stdev, to be continued
            avg_terminal_branch_len = sum_terminal_branch_len / len(terminals)
            tree_value = len(internals) / len(terminals)
            clean()
    sys.setrecursionlimit(old_max_recursion)
    return (gap_ratio, resolution, entropy, pi, tree_value,
                avg_terminal_branch_len)


def evaluate_main(arg_str):
    """
    Evaluate variance of alignments.
    Args:
        arg_str:

    Returns:
    """
    if arg_str is None:
        arg = parse_args()
    else:
        arg = parse_args(arg_str.split(' '))