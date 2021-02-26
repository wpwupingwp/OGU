#!/usr/bin/python3

import argparse
import logging
import sys
import numpy as np

from collections import namedtuple
from io import StringIO
from os import devnull, cpu_count
from pathlib import Path
from subprocess import run

from Bio import Phylo, SeqIO
from matplotlib import use as mpl_use
mpl_use('Agg')
from matplotlib import pyplot as plt
from matplotlib import rcParams
params = {'axes.labelsize': 12, 'axes.linewidth': 1.5, 'axes.titlesize': 20,
          'font.size': 12, 'lines.linewidth': 1.5,
          'legend.fontsize': 10, 'legend.handlelength': 2}
rcParams.update(params)
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


def parse_args():
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=evaluate_main.__doc__)
    arg.add_argument('-fasta', nargs='*', help='unaligned fasta files')
    arg.add_argument('-fasta_folder', default=None, help='folder of fasta files')
    arg.add_argument('-aln', nargs='*', help='aligned files')
    arg.add_argument('-out', help='output folder')
    arg.add_argument('-skip_primer', action='store_true',
                     help='skip primer designing')
    options = arg.add_argument_group('Options')
    options.add_argument('-ig', '-ignore_gap', dest='ignore_gap',
                         action='store_true',
                         help='ignore gaps in alignments')
    options.add_argument('-iab', '-ignore_ambiguous',
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
    return arg.parse_known_args()


def init_arg(arg):
    arg = utils.init_out(arg)
    if arg.out is None:
        return None
    if all([arg.fasta, arg.fasta_folder]):
        log.info('Do not recommend to  use "-fasta" and "-fasta_folder" '
                 'at same time!')
    if arg.fasta is not None:
        arg.fasta = [Path(i).absolute() for i in arg.fasta]
    else:
        arg.fasta = []
    if arg.fasta_folder is not None:
        for i in Path(arg.fasta_folder).glob('*'):
            arg.fasta.append(i.absolute())
    log.info('Add fasta from gb2fasta module if possible.')
    for i in Path(arg._unique).glob('*'):
        if i.name != 'Unknown.fasta':
            arg.fasta.append(i.absolute())
    if arg.fasta:
        for i in arg.fasta:
            if not i.exists() or not i.is_file():
                log.error(f'{i} does not exist or is not a valid file.')
                return None
    if arg.aln is not None:
        arg.aln = [Path(i).absolute() for i in arg.aln]
        for i in arg.aln:
            if not i.exists() or not i.is_file():
                log.error(f'{i} does not exist or is not a valid file.')
                return None
    else:
        arg.aln = []
    if not any([arg.fasta, arg.aln, arg.fasta_folder]):
        log.warning('Empty input.')
        return None
    if arg.quick:
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
    if not files:
        return aligned, unaligned
    _, mafft = utils.get_mafft()
    if not _:
        log.error('Cannot run MAFFT.')
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
        log.error(f'Invalid alignment file {aln_fasta}')
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
    n_columns = alignment.shape[1]
    n_gap_columns = gap_columns.shape[1]
    log.info(f'{n_columns} columns, {n_gap_columns} have gaps.')
    if n_gap_columns/n_columns > 0.5:
        log.warning('Too much columns with gaps.')
    return no_gap_columns, gap_columns


def old_remove_gap(aln_fasta: Path, new_file: Path) -> Path:
    # to be removed
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
    SeqIO.convert(no_gap, 'fasta', new_file, 'fasta')
    no_gap.close()
    return new_file


def gc_ratio(alignment: np.array, ignore_ambiguous=True) -> (
        float, np.array):
    """
    Get GC ratio of total alignment and each sequence.
    Args:
        alignment: np.array
        ignore_ambiguous: count ambiguous bases or not
    Returns:
        total_gc: gc value
        gc_array: np.array of gc value for each row(sequence)
    """
    def get_gc_ratio(row: np.array, ignore: bool):
        # if True, do not count ambiguous bases as GC, but count them in total
        # BDVH = 1/3 GC
        # MSYKRS = 1/2 GC
        gc = 0
        count_ = np.unique(row, return_counts=True)
        count = dict(zip(count_[0], count_[1]))
        for base in (b'G', b'C'):
            gc += count.get(base, 0)
        if not ignore:
            gc += count.get(b'N', 0) / 4
            for ambiguous_base in (b'B', b'D', b'V', b'H'):
                gc += count.get(ambiguous_base, 0) / 3
            for ambiguous_base in (b'M', b'S', b'Y', b'K', b'R', b'S'):
                gc += count.get(ambiguous_base, 0) / 2
        return gc / (rows*columns-count.get(b'-', 0))

    rows, columns = alignment.shape
    total_gc = get_gc_ratio(alignment, ignore_ambiguous)
    gc_array = np.fromiter(
        [get_gc_ratio(row, ignore_ambiguous) for row in alignment],
        dtype=np.float)
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
    # pi should > 0, use max()
    return max(0, pi)


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
        pd_stem: only calculate stem branch
        pd_stem_sd: standard deviation
        pd_terminal: only calculate terminal
        pd_terminal_sd: standard deviation
        tree_res: tree resolution
    """
    pd = 0.0
    pd_stem = 0.0
    pd_stem_sd = 0.0
    pd_terminal = 0.0
    pd_terminal_sd = 0.0
    tree_res = 0.0
    rows, columns = alignment.shape
    if rows < 4:
        log.debug('Too few sequences.')
        return pd, pd_stem, pd_stem_sd, pd_terminal, pd_terminal_sd, tree_res
    old_max_recursion = sys.getrecursionlimit()
    sys.setrecursionlimit(max(rows+10, old_max_recursion))
    aln_file = tmp / f'{columns}.tmp'
    array_to_fasta(alignment, aln_file)
    _, iqtree = utils.get_iqtree()
    if not _:
        log.critical('Cannot find iqtree.')
        return pd, pd_stem, pd_stem_sd, pd_terminal, pd_terminal_sd, tree_res
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
            pd_stem = sum([i.branch_length for i in internals])
            pd_stem_sd = np.std([i.branch_length for i in internals])
            pd_terminal = sum([i.branch_length for i in terminals])
            pd_terminal_sd = np.std([i.branch_length for i in terminals])
            # avg_terminal_branch_len = pd_terminal / rows
            # may be zero
            tree_res = len(internals) / max(1, len(terminals))
        except Exception:
            log.info('Bad phylogenetic tree.')
    utils.clean_tmp(aln_file)
    sys.setrecursionlimit(old_max_recursion)
    return pd, pd_stem, pd_stem_sd, pd_terminal, pd_terminal_sd, tree_res


class Variance(namedtuple('Variance',
                          ['Samples', 'Length', 'Gap_Ratio',
                           'Observed_Res', 'Entropy', 'Pi',
                           'PD', 'PD_stem', 'PD_stem_SD',
                           'PD_terminal', 'PD_terminal_SD',
                           'Tree_Res', 'Total_GC'],
                          defaults=[0 for i in range(13)])):
    """
    For get_resolution()
    samples: rows
    length: columns
    """
    __slots__ = ()

    def __str__(self):
        return ('{Samples},{Length},{Gap_Ratio:.4%},'
                '{Observed_Res:.4%},{Entropy:.8f},{Pi:.8f},'
                '{PD:.8f},{PD_stem:.8f},{PD_stem_SD:.8f},'
                '{PD_terminal:.8f},{PD_terminal_SD:.8f},'
                '{Tree_Res:.4%},{Total_GC:.4%}'.format(**self._asdict()))


def get_resolution(alignment: np.array, tmp: Path,
                   ignore_ambiguous=True) -> tuple:
    """
    Given alignment (2d numpy array), location of fragment(start and end, int,
    start from zero, exclude end),
    return gap ratio, resolution, entropy, Pi, tree value and average terminal
    branch length.
    """
    gc_array = np.array([0])
    rows, columns = alignment.shape
    total = rows * columns
    # index error
    if columns == 0:
        return Variance(), gc_array
    gap_ratio = len(alignment[alignment == b'-']) / total
    item, count = np.unique(alignment, return_counts=True, axis=0)
    observed_res = len(count) / rows
    # normalized entropy
    entropy = normalized_entropy(count, rows)
    pi = nucleotide_diversity(alignment)
    (pd, pd_stem, pd_stem_sd, pd_terminal, pd_terminal_sd,
     tree_res) = phylogenetic_diversity(alignment, tmp)
    total_gc, gc_array = gc_ratio(alignment, ignore_ambiguous)
    variance = Variance(rows, columns, gap_ratio, observed_res, entropy, pi,
                        pd, pd_stem, pd_stem_sd, pd_terminal, pd_terminal_sd,
                        tree_res, total_gc)
    return variance, gc_array


def output_sliding(sliding: list, name: str, out: Path,
                   size: int, step: int) -> (Path, Path):
    if len(sliding) == 0:
        log.warning('Empty sliding-window result.')
        return Path(), Path()
    out_csv = out / (name+'.csv')
    head = 'Start,End,' + ','.join(Variance._fields) + '\n'
    handle = open(out_csv, 'w', encoding='utf-8')
    handle.write(head)
    index = []
    # for human reading
    start = 1
    for variance in sliding:
        line = f'{start},{start+size},{variance}\n'
        index.append(start)
        handle.write(line)
        start += step
    handle.close()
    # draw
    out_pdf = out / (name+'.pdf')
    plt.style.use('seaborn-colorblind')
    # how to find optimized size?
    fig, ax1 = plt.subplots(figsize=(15 + len(sliding) // 5000, 10))
    ax1.yaxis.set_ticks(np.linspace(0, 1, num=11))
    ax1.set_title(f'Sliding window results of {name} '
                  f'(sample={sliding[0].Samples}, '
                  f'size={size} bp, step={step} bp)', pad=30)
    ax1.set_xlabel('Bases')
    ax1.set_ylabel('Gap Ratio, GC Ratio, Resolution & Shannon Index')
    ax1.plot(index, [i.Gap_Ratio for i in sliding], label='Gap Ratio',
             alpha=0.8)
    ax1.plot(index, [i.Observed_Res for i in sliding],
             label='Observed Resolution', alpha=0.8)
    ax1.plot(index, [i.Entropy for i in sliding],
             label='Shannon Equitability Index', alpha=0.8)
    ax1.plot(index, [i.Tree_Res for i in sliding], label='Tree Resolution',
             alpha=0.8)
    ax1.plot(index, [i.Total_GC for i in sliding], label='GC Ratio',
             alpha=0.8)
    ax1.legend(loc='lower left')
    # different ytick
    ax2 = ax1.twinx()
    # ax2..yaxis.set_ticks(np.linspace(0, max_range, 21))
    ax2.set_ylabel(r'$\pi$ & Phylogenetic Diversity', rotation=-90, labelpad=20)
    ax2.plot(index, [i.PD for i in sliding], linestyle='--', label='PD',
             alpha=0.8)
    ax2.plot(index, [i.PD_stem for i in sliding], linestyle='--',
             label='PD_stem', alpha=0.8)
    ax2.plot(index, [i.PD_terminal for i in sliding], linestyle='--',
             label='PD_terminal', alpha=0.8)
    # k--
    ax2.plot(index, [i.Pi for i in sliding], '--', label=r'$\pi$', alpha=0.8)
    ax2.legend(loc='upper right')
    plt.savefig(out_pdf)
    plt.close()
    return out_csv, out_pdf


def evaluate(aln: Path, arg) -> tuple:
    """
    Wrapper
    Args:
        aln: alignment file
        arg: args
    Returns:
        summary: namedtuple
        gc_array: for draw
        sliding: list of Variance
    """
    sliding = []
    name, alignment = fasta_to_array(aln)
    if name is None:
        log.info(f'Invalid fasta file {aln}.')
        return None, None, None
    if arg.ignore_gap:
        no_gap_alignment, gap_alignment = remove_gap(alignment)
    else:
        no_gap_alignment = alignment
        gap_alignment = np.array([[]])
    rows, columns = alignment.shape
    log.info(f'Evaluate {aln}')
    summary, gc_array = get_resolution(no_gap_alignment, arg._tmp,
                                       arg.ignore_ambiguous_base)
    if arg.quick:
        pass
    else:
        # sliding window
        for i in range(0, columns, arg.step):
            # view, not copy
            subalign = alignment[:, i:i+arg.size]
            variance, sub_gc_array = get_resolution(subalign, arg._tmp,
                                                    arg.ignore_ambiguous_base)
            sliding.append(variance)
    return summary, gc_array, sliding


def evaluate_main():
    """
    Evaluate variance of alignments.
    """
    log.info('Running evaluate module...')
    arg, other_args2 = parse_args()
    arg = init_arg(arg)
    if arg is None:
        log.info('Quit evaluate module.')
        return None, other_args2
    utils.add_file_log(arg)
    aligned, unaligned = align(arg.fasta, arg._align)
    aligned.extend(arg.aln)
    evaluation_result = arg.out / 'Evaluation.csv'
    csv_head = 'Loci,' + ','.join(Variance._fields) + '\n'
    with open(evaluation_result, 'w', encoding='utf-8') as out_csv:
        out_csv.write(csv_head)
    for aln in aligned:
        summary, gc_array, sliding = evaluate(aln, arg)
        if summary is None:
            continue
        log.info(f'\tGap ratio:                 {summary.Gap_Ratio:.8f}')
        log.info(f'\tObserved resolution:       {summary.Observed_Res:.8f}')
        log.info(f'\tNormalized Shannon Index:  {summary.Entropy:.8f}')
        log.info(f'\tPi:                        {summary.Pi:.8f}')
        log.info(f'\tPhylogenetic diversity:    {summary.PD:.8f}')
        log.info(f'\tTree resolution:           {summary.Tree_Res:.8f}')
        with open(evaluation_result, 'a', encoding='utf-8') as out:
            out.write(aln.stem+','+str(summary)+'\n')
        if not arg.quick:
            output_sliding(sliding, aln.stem, arg._evaluate, arg.size, arg.step)
    log.info(f'Evaluation results could be found in {evaluation_result}')
    log.info('Evaluate module finished.')
    # for i in aligned:
    #     utils.move(i, arg._align/(i.name), copy=True)
    return arg, other_args2


if __name__ == '__main__':
    evaluate_main()