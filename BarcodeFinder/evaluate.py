#!/usr/bin/python3

import argparse
import logging
import sys
import numpy as np

from io import StringIO
from os import devnull, cpu_count
from pathlib import Path
from subprocess import run

from Bio import Phylo

from BarcodeFinder import utils

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
    options.add_argument('-ig', '-ignore_gap', dest='ignore_gap',
                         action='store_true',
                         help='ignore gaps in alignments')
    options.add_argument('-iab', '-ignore_ambiguous_base',
                         dest='ignore_ambiguous_base',
                         action='store_true',
                         help='ignore ambiguous bases like "M" or "N"')
    sliding_window = arg.add_argument_group('Sliding-window')
    sliding_window.add_argument('-quick', action='store_true',
                                help='skip sliding-window analysis')
    sliding_window.add_argument('-size', type=int, default=500,
                                help='window size')
    sliding_window.add_argument('-step', type=int, default=50,
                                help='step length')
    if arg_list is None:
        return arg.parse_args()
    else:
        return arg.parse_args(arg_list)


def init_arg(arg):
    if arg.fasta is None and arg.aln is None:
        log.error('Empty input.')
        return None
    if arg.fasta is not None:
        arg.fasta = [Path(i).absolute() for i in arg.fasta]
    if arg.aln is not None:
        arg.aln = [Path(i).absolute() for i in arg.aln]
    arg.out = utils.init_out(arg)
    if arg.fast:
        log.info('The "-quick" mode was opened. '
                 'Skip sliding-window analysis')
    return arg


def align(files: list, folder: Path) -> (list, list):
    """
    Align sequences with mafft.
    Args:
        files(list): fasta files
        folder(path): folder for output
    Returns:
        aligned(list): aligned fasta
        unaligned(list): unaligned files

    """
    log.info('Align sequences.')
    aligned = list()
    unaligned = list()
    _, mafft = utils.get_mafft()
    if not _:
        log.error('Cannot run mafft.')
        return aligned, unaligned
    # get available CPU cores
    cores = max(1, cpu_count() - 1)
    for fasta in files:
        log.info('Aligning {}.'.format(fasta))
        out = folder / fasta.with_suffix('.aln').name
        with open(devnull, 'w', encoding='utf-8') as f:
            # if computer is good enough, "--genafpair" is recommended
            # where is mafft?
            _ = (f'{mafft} --auto --thread {cores} --reorder --quiet '
                 f'--adjustdirection {fasta} > {out}')
            m = run(_, shell=True, stdout=f, stderr=f)
        if m.returncode == 0:
            aligned.append(out)
        else:
            unaligned.append(fasta)
            try:
                out.unlink()
            except Exception:
                pass
    log.info('Alignment finished.')
    for i in Path().cwd().glob('_order*'):
        i.unlink()
    log.info(f'Total {len(files)} files')
    log.info(f'Aligned {len(aligned)}')
    log.info(f'Unaligned {len(unaligned)}')
    return aligned, unaligned


def remove_gap(alignment: np.array) -> (np.array, np.array):
    """
    Split alignment into with_gap and without_gap.
    Args:
        alignment: raw array
    Returns:
        no_gap_columns: without gap
        gap_columns: columns having gaps
    """
    gap = b'-'
    # axis 0 for column
    have_gap = np.any(alignment==gap, axis=0)
    gap_columns = alignment[:, have_gap]
    no_gap_columns = alignment[:, ~have_gap]
    n_columns = alignment.shape()[1]
    n_gap_columns = gap_columns.shape()[1]
    log.info(f'{n_columns} columns, {n_gap_columns} have gaps.')
    if n_gap_columns/n_columns > 0.5:
        log.warning('Too much columns with gaps.')
    return no_gap_columns, gap_columns


def old_remove_gap(aln_fasta: Path, new_file: Path) -> Path:
    """
    old function, for BLAST
    Args:
        aln_fasta: fasta with gap
        new_file: fasta without gap
    Returns:
        new_file: fasta without gap
    """
    no_gap = StringIO()
    with open(aln_fasta, 'r', encoding='utf-8') as raw:
        for line in raw:
            no_gap.write(line.replace('-', ''))
    # try to avoid makeblastdb error
    no_gap.seek(0)
    from Bio import SeqIO
    SeqIO.convert(no_gap, 'fasta', new_file, 'fasta')
    no_gap.close()
    return new_file


def array_to_fasta(alignment: np.array, filename: Path) -> Path:
    """
    Convert np.array to fasta.
    Use index number as sequence id.
    Args:
        alignment
        filename
    Returns:
        filename
    """
    with open(filename, 'wb') as aln:
        for index, row in enumerate(alignment):
            aln.write(b'>'+str(index).encode('utf-8')+b'\n')
            aln.write(b''.join(row)+b'\n')
    return filename


def fasta_to_array(aln_fasta: Path) -> (np.array, np.array):
    """
    Given fasta format alignment filename, return a numpy array for sequence:
    Faster and use smaller mem.
    Ensure all bases are capital.
    Args:
        aln_fasta(Path): aligned fasta file
    Returns:
        name(np.array): name array
        sequence(np.array): sequence array
    """
    data = []
    record = ['id', 'sequence']
    with open(aln_fasta, 'r', encoding='utf-8') as raw:
        for line in raw:
            if line.startswith('>'):
                data.append([record[0], ''.join(record[1:])])
                # remove ">" and CRLF
                name = line[1:].strip()
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
        return None, None
    # Convert List to numpy array.
    # order 'F' is a bit faster than 'C'
    # new = np.hstack((name, seq)) -> is slower
    name_array = np.array([[i[0]] for i in data], dtype=np.bytes_)
    # fromiter is faster than from list
    # S1: bytes
    sequence_array = np.array(
        [np.fromiter(i[1], dtype=np.dtype('S1')) for i in data],
        order='F')
    if name_array is None:
        log.error('Bad fasta file {}.'.format(aln_fasta))
    return name_array, sequence_array


def gc_ratio(alignment: np.array, ignore_ambiguous_base=True) -> (
        float, np.array):
    def get_gc_ratio(row: np.array, ignore: True):
        # if True, do not count ambiguous bases as GC, but count them in total
        # BDVH = 1/3 GC
        # MSYKRS = 1/2 GC
        gc = 0
        count = np.unique(row, return_counts=True)
        for base in b'GC':
            gc += count[base]
        if not ignore_ambiguous_base:
            for ambiguous_base in b'BDVH':
                gc += count[ambiguous_base] / 3
            for ambiguous_base in b'MSYKRS':
                gc += count[ambiguous_base] / 2
        return gc / (rows*columns-counts[b'-'])

    rows, columns = alignment.shape()
    total_gc = get_gc_ratio(alignment, ignore_ambiguous_base)
    gc_array = np.fromiter([
        get_gc_ratio(row, ignore_ambiguous_base) for row in alignment])
    return total_gc, gc_array


def normalized_entropy(count: np.array, rows: int) -> float:
    """
    Calculate normalized entropy.
    Args:
        count: np.unique(axis=0)
        rows: rows number
    Returns:
        entropy(float): normalized entropy
    """
    max_h = np.log2(rows)
    entropy = 0
    for j in count:
        p_j = j / rows
        log2_p_j = np.log2(p_j)
        entropy += log2_p_j * p_j
    entropy = -1 * entropy / max_h
    # entropy should > 0
    return max(0, entropy)


def nucleotide_diversity(alignment: np.array) -> float:
    """
    Nucleotide diversity (pi)
    Args:
        alignment: np.array
    Returns:
        pi: float
    """
    rows, columns = alignment.shape
    m = columns
    n = rows
    sum_d_ij = 0
    for i in range(n):
        d_ij = np.sum(alignment[i] != alignment[(i + 1):])
        sum_d_ij += d_ij
    pi = (2 / (n * (n - 1)) * sum_d_ij) / m
    # pi should > 0
    return min(0, pi)


def phylogenetic_diversity(alignment: np.array, tmp: Path) -> (float, float,
                                                               float, float):
    """
    Calculate the phylogenetic diversity.
    Use HKY model for saving time.
    Args:
        alignment: np.array
        tmp: tmp folder
    Returns:
        pd: phylogenetic diversity
        pd_terminal: only calculate terminal
        pd_stem: only calculate stem branch
        tree_res: tree resolution
    """
    pd = 0.0
    pd_terminal = 0.0
    pd_stem = 0.0
    tree_res = 0.0
    rows, columns = alignment.shape()
    if rows < 4:
        log.debug('Too few sequences.')
        return pd, pd_stem, pd_terminal, tree_res
    old_max_recursion = sys.getrecursionlimit()
    sys.setrecursionlimit(max(rows+10, old_max_recursion))
    aln_file = tmp / f'{columns}.tmp'
    array_to_fasta(alignment, aln_file)
    _, iqtree = utils.get_iqtree()
    if not _:
        log.critical('Cannot find iqtree.')
        return pd, pd_terminal, pd_stem, tree_res
    with open(devnull, 'w', encoding='utf-8') as out:
        run_ = run(f'{iqtree} -s {aln_file} -m HKY -fast -czb -redo',
                   stdout=out, stderr=out, shell=True)
    # just return 0 if there is error
    if run_.returncode != 0:
        log.debug('Too much gap in the alignment.')
    else:
        tree = Phylo.read(str(aln_file)+'.treefile', 'newick')
        # skip the first empty node
        try:
            pd = tree.total_branch_length()
            internals = tree.get_nonterminals()[1:]
            terminals = tree.get_terminals()
            pd_terminal = sum([i.branch_length for i in terminals])
            pd_stem = sum([i.branch_length for i in internals])
            # miss stdev, to be continued
            avg_terminal_branch_len = pd_terminal / rows
            # may be zero
            tree_res = len(internals) / max(1, len(terminals))
        except Exception:
            log.info('Bad phylogenetic tree.')
    utils.clean_tmp(aln_file)
    sys.setrecursionlimit(old_max_recursion)
    return pd, pd_terminal, pd_stem, tree_res


def get_resolution(alignment: np.array, start: int, end: int,
                   tmp: Path, quick=False) -> tuple:
    """
    Given alignment (2d numpy array), location of fragment(start and end, int,
    start from zero, exclude end),
    return gap ratio, resolution, entropy, Pi, tree value and average terminal
    branch length.
    """
    subalign = alignment[:, start:end]
    rows, columns = subalign.shape
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
    entropy = normalized_entropy(count, rows)
    pi = nucleotide_diversity(subalign)
    # Nucleotide diversity (pi)
    # tree value
    pd, pd_stem, pd_terminal, tree_res = phylogenetic_diversity(subalign, tmp)
    return (gap_ratio, resolution, entropy, pi, pd, pd_stem, pd_terminal,
            tree_res)


def gap_analyze(gap_alignment: np.array):
    pass


def evaluate(aln: Path, result: Path, arg) -> bool:
    name, alignment = fasta_to_array(aln)
    if name is None:
        log.info('Invalid fasta file {}.'.format(aln))
        return False
    if arg.ignore_gap:
        no_gap_alignment, gap_alignment = remove_gap(alignment)
    else:
        no_gap_alignment = alignment
        gap_alignment = np.array([[]])
    rows, columns = alignment.shape
    log.info(f'Evaluate {aln}')
    (gap_ratio, observed_res, entropy, pi, tree_res,
     branch_len) = get_resolution(no_gap_alignment, 0, columns, arg.tmp)
    log.info(f'\tGap ratio:\t{gap_ratio:.8f}')
    log.info(f'\tObserved resolution:\t{observed_res:.8f}')
    log.info(f'\tNormalized Shannon Index:\t{entropy:.8f}')
    log.info(f'\tPi:\t{pi:.8f}')
    log.info(f'\tTree resolution:\t{tree_res:.8f}')
    log.info(f'\tAverage terminal branch length:\t{branch_len:.8f}')
    with open(result, 'a', encoding='utf-8') as out:
        out.write('{},{},{},{:.4%},{:.8f},{:.8f},{:.8f},{:.8f},{:.8f}'
                  '\n'.format(aln.stem, rows, columns, gap_ratio,
                              observed_res, tree_res, entropy,
                              branch_len, pi))
    return True


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
    if arg is None:
        log.info('Quit.')
        return None
    aligned, unaligned = align(arg.fasta, arg._align)
    aligned.extend(arg.aln)
    evaluation_result = arg.out / 'Evaluation.csv'
    with open(evaluation_result, 'w', encoding='utf-8') as out_csv:
       out_csv.write('Loci,Samples,Length,GapRatio,ObservedResolution,'
                     'TreeResolution,ShannonIndex,AvgTerminalBranchLen,'
                     'Pi\n')
    for aln in aligned:
        evaluate(aln, evaluation_result, arg)