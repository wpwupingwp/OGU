#!/usr/bin/python3

import argparse
import os
import re

import numpy as np
import primer3

from collections import defaultdict
from math import log2
from multiprocessing import cpu_count
from timeit import default_timer as timer
from subprocess import run
from typing import List, Dict, Any

from Bio import SearchIO, SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline as nb
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

import matplotlib
from matplotlib import pyplot as plt
from matplotlib import ticker as mtick
matplotlib.rcParams['lines.linewidth'] = 1.5
matplotlib.rcParams['axes.linewidth'] = 1.5
matplotlib.rcParams['axes.labelsize'] = 16
matplotlib.rcParams['axes.titlesize'] = 25
matplotlib.rcParams['font.size'] = 16
matplotlib.rcParams['axes.facecolor'] = '#666666'
matplotlib.rcParams['figure.figsize'] = 16, 9


class PrimerWithInfo(SeqRecord):
    def __init__(self, seq='', quality='', index=0, start=0, tm=0,
                 coverage=0, sum_bitscore=0, avg_mid_loc=0, detail=0):
        super().__init__(seq=seq)

        self.sequence = self.seq
        self.quality = self.letter_annotations['solexa_quality'] = quality
        self.index = self.annotations['index'] = index
        self.start = self.annotations['start'] = start
        self.tm = self.annotations['tm'] = tm
        self.coverage = self.annotations['coverage'] = coverage
        self.sum_bitscore = self.annotations['sum_bitscore'] = sum_bitscore
        self.avg_mid_loc = self.annotations['avg_mid_loc'] = avg_mid_loc
        self.detail = self.annotations['detail'] = detail
        self.end = self.annotations['end'] = start + self.__len__() - 1
        self.update_id()

    def update_id(self):
        self.id = ('{}-Start({})-End({})-Tm({:.2f})-Coverage({:.2%})-'
                   'SumBitScore({})-AvgMidLocation({:.0f})'.format(
                      self.index, self.start, self.end, self.tm, self.coverage,
                      self.sum_bitscore, self.avg_mid_loc))

    def __getitem__(self, i):
        if isinstance(i, int):
            i = slice(i, i+1)
        if isinstance(i, slice):
            if self.seq is None:
                raise ValueError('Empty sequence')
            answer = PrimerWithInfo(seq=self.seq[i], quality=self.quality[i])
            answer.annotations = dict(self.annotations.items())
            return answer

#@profile
def prepare(fasta):
    """
    Given fasta format alignment filename, return a numpy array for sequence:
    Generate fasta file without gap for makeblastdb, return file name.
    """
    no_gap = fasta + '.no_gap'
    data: List[List[str, str]] = []
    record = ['id', 'sequence']
    with open(fasta, 'r') as raw, open(no_gap, 'w') as out:
        for line in raw:
            out.write(line.replace('-', ''))
            if line.startswith('>'):
                data.append([record[0], ''.join(record[1:])])
                # remove ">" and CRLF
                name = line.strip('>\r\n')
                record = [name, '']
            else:
                record.append(line.strip().upper())
        # add last sequence
        data.append([record[0], ''.join(record[1:])])
        # generate no-gap fasta
    # skip head['id', 'seq']
    data = data[1:]
    # check sequence length
    length_check = [len(i[1]) for i in data]
    assert len(set(length_check)) == 1, (
        'Alignment does not have uniform width, please check again !')

    # Convert List to numpy array.
    # order 'F' is a bit faster than 'C'
    # new = np.hstack((name, seq)) -> is slower
    name = np.array([[i[0]] for i in data], dtype=np.bytes_)
    sequence = np.array([list(i[1]) for i in data], dtype=np.bytes_, order='F')
    return name, sequence, no_gap


#@profile
def count_base(alignment, rows, columns):
    """
    Given alignment numpy array, count cumulative frequency of base in each
    column (consider ambiguous base and "N", "-" and "?", otherwise omit).
    Return List[List[float, float, float, float, float, float, float]] for
    [A, T, C, G, N, GAP, OTHER].
    """
    frequency: List[List[float, float, float, float, float, float, float]] = []
    for index in range(columns):
        column = alignment[:, [index]]
        base, counts = np.unique(column, return_counts=True)
        count_dict = {b'A': 0, b'C': 0, b'G': 0, b'T': 0, b'M': 0, b'R': 0,
                      b'W': 0, b'S': 0, b'Y': 0, b'K': 0, b'V': 0, b'H': 0,
                      b'D': 0, b'B': 0, b'X': 0, b'N': 0, b'-': 0, b'?': 0}
        count_dict.update(dict(zip(base, counts)))
        a = (count_dict[b'A'] +
             (count_dict[b'D']+count_dict[b'H']+count_dict[b'V'])/3 +
             (count_dict[b'M']+count_dict[b'R'] + count_dict[b'W'])/2)
        t = (count_dict[b'T'] +
             (count_dict[b'B']+count_dict[b'H']+count_dict[b'D'])/3 +
             (count_dict[b'K']+count_dict[b'W'] + count_dict[b'Y'])/2)
        c = (count_dict[b'C'] +
             (count_dict[b'B']+count_dict[b'H']+count_dict[b'V'])/3 +
             (count_dict[b'M']+count_dict[b'S'] + count_dict[b'Y'])/2)
        g = (count_dict[b'G'] +
             (count_dict[b'B']+count_dict[b'D']+count_dict[b'V'])/3 +
             (count_dict[b'K']+count_dict[b'R'] + count_dict[b'S'])/2)
        gap = count_dict[b'-']
        n = count_dict[b'N'] + count_dict[b'X'] + count_dict[b'?']
        other = rows - a - t - c - g - gap - n
        frequency.append([a, t, c, g, n, gap, other])
    return frequency


#@profile
def get_quality_string(data: List[float], rows: int):
    # use fastq-illumina format
    max_q = 62
    factor = max_q/rows
    # use min to avoid KeyError
    quality_value = [min(max_q, int(i*factor))-1 for i in data]
    return quality_value 


#@profile
def generate_consensus(base_cumulative_frequency, cutoff, gap_cutoff,
                       rows, columns, output):
    """
    Given base count info, return List[index, base, quality]
    and List[List[str, str, str, PrimerInfo]] for writing conesensus.
    """
    def get_ambiguous_dict():
        from Bio.Data.IUPACData import ambiguous_dna_values
        data = ambiguous_dna_values
        data = dict(zip(data.values(), data.keys()))
        # 2:{'AC': 'M',}
        data_with_len = defaultdict(lambda: dict())
        for key in data:
            data_with_len[len(key)][key] = data[key]
        return data_with_len

    ambiguous_dict = get_ambiguous_dict()

    most: List[List[int, str, float]] = []
    gap_cutoff = rows * gap_cutoff
    cutoff = rows * cutoff

    for location, column in enumerate(base_cumulative_frequency, 1):
        finish = False
        # "*" for others
        value = dict(zip(list('ATCGN-*'), column))

        base = 'N'
        if value['N'] >= gap_cutoff:
            count = value['N']
            most.append([location, base, count])
            continue
        sum_gap = sum([value['-'], value['*']])
        if sum_gap >= gap_cutoff:
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
                    if count >= cutoff:
                        base = ambiguous_dict[length][key]
                        finish = True
                        most.append([location, base, count])
    quality = [i[2] for i in most]
    consensus = PrimerWithInfo(start=1, seq=''.join([i[1] for i in most]),
                               quality=get_quality_string(quality, rows))
    print(len(consensus), len(quality))
    SeqIO.write(consensus, output, 'fastq')
    return consensus


#@profile
def find_continuous(consensus, good_region, min_len):
    """
    Given PrimerWithInfo, good_region: List[bool], min_len
    Return consensus with features.
    """
    skip = ('N', '-')
    start = 0
    for index, base in enumerate(consensus.sequence[-min_len:]):
        if not good_region[index] or base in skip:
            if (index-start) >= min_len:
                consensus.features.append(SeqFeature( FeatureLocation(
                    start, index), strand=1))
            start = index
    for i in consensus.features:
        print(i)
    return consensus


#@profile
def find_primer(continuous, rows, min_len, max_len, ambiguous_base_n):
    """
    Find suitable primer in given List[List[int, str, float]]
    return List[List[str, str, PrimerInfo]]
    """
    poly = re.compile(r'([ATCG])\1\1\1\1')
    ambiguous_base = re.compile(r'[^ATCG]')
    tandem = re.compile(r'([ATCG]{2})\1\1\1\1')

    def is_good_primer(primer):
        # ref1. http://www.premierbiosoft.com/tech_notes/PCR_Primer_Design.html
        seq = ''.join(primer)
        if re.search(poly, seq) is not None:
            return False, 0, 'Poly(NNNNN) structure found'
        if re.search(tandem, seq) is not None:
            return False, 0, 'Tandom(NN*5) exist'
            # no more 3 ambiguous base
        if len(re.findall(ambiguous_base, seq)) >= ambiguous_base_n:
            return False, 0, 'More than 3 ambiguous base'

# primer3.setGlobals seems have no effect on calcTm, so I have to replace all
# ambiguous base to A to get an approximate value. Othervise calcTm() will
# generate -99999 if there is ambiguous base.
        pure_seq = re.sub(ambiguous_base, 'A', seq)
        tm = primer3.calcTm(pure_seq)
        hairpin_tm = primer3.calcHairpinTm(pure_seq)
        homodimer_tm = primer3.calcHomodimerTm(pure_seq)
        if max(tm, hairpin_tm, homodimer_tm) != tm:
            return False, 0, 'Hairpin or homodimer found'
        return True, tm, 'Ok'

    primers: List[PrimerWithInfo] = list()
    continuous = [i for i in continuous if len(i) >= min_len]
    n = 1
    for fragment in continuous:
        print(fragment)
        len_fragment = len(fragment)
        for begin in range(len_fragment-max_len):
            for p_len in range(min_len, max_len):
                seqs = fragment[begin:(begin+p_len)]
                good_primer, tm, detail = is_good_primer(seqs)
                if good_primer:
                    start = seqs[0][0]
                    sequence = ''.join([i[1] for i in seqs])
                    quality = [i[2] for i in seqs]
                    primers.append(PrimerWithInfo(
                        index=n, start=start, tm=tm, sequence=sequence,
                        quality=get_quality_string(quality, rows)))
                    n += 1
                else:
                    continue
    return primers


#@profile
def count_and_draw(alignment, consensus, arg):
    """
    Given alignment(numpy array), return unique sequence count List[float].
    Calculate Shannon Index based on
    http://www.tiem.utk.edu/~gross/bioed/bealsmodules/shannonDI.html
    return List[float]
    All calculation excludes primer sequence.
    """
    rows, columns = alignment.shape
    min_primer = arg.min_primer
    max_primer = arg.max_primer
    min_product = arg.min_product
    max_product = arg.max_product
    window = arg.window
    out = arg.out
    # Different count
    count_min_len: List[List[float]] = list()
    count_max_len: List[List[float]] = list()
    shannon_index1: List[float] = list()
    shannon_index2: List[float] = list()
    max_shannon_index = -1*((1/rows)*log2(1/rows)*rows)
    index: List[int] = list()
    min_plus = min_product - max_primer * 2
    max_plus = max_product - min_primer * 2
    for i in range(0, columns-min_product, window):
        # skip gap
        if consensus.sequence[i] in ('-', 'N'):
            continue
        # exclude primer sequence
        split_1 = alignment[:, i:(i+min_plus)]
        split_2 = alignment[:, i:(i+max_plus)]
        # uniqe array, count line*times
        _, count1 = np.unique(split_1, return_counts=True, axis=0)
        _, count2 = np.unique(split_2, return_counts=True, axis=0)
        count_min_len.append(len(count1))
        count_max_len.append(len(count2))
        # Shannon Index
        h = 0
        for j in count1:
            p_j = j / rows
            log2_p_j = log2(p_j)
            h += log2_p_j * p_j
        shannon_index1.append(-1*h)
        h = 0
        for j in count2:
            p_j = j / rows
            log2_p_j = log2(p_j)
            h += log2_p_j * p_j
        shannon_index2.append(-1*h)
        index.append(i)

    # convert value to (0,1)
    count_min_len = [i/rows for i in count_min_len]
    count_max_len = [i/rows for i in count_max_len]
    # plt.style.use('ggplot')
    fig, ax1 = plt.subplots()
    plt.title('Shannon Index & Resolution({}-{}bp, window={})'.format(
        min_product, max_product, window))
    plt.xlabel('Base')
    plt.xticks(range(0, columns, int(columns/10)))
    # c=List for different color, s=size for different size
    ax1.scatter(index, shannon_index1, c=shannon_index1,
                cmap='GnBu', alpha=0.8, s=10, label='{}bp'.format(min_product))
    ax1.scatter(index, shannon_index2, c=shannon_index2,
                cmap='OrRd', alpha=0.8, s=10, label='{}bp'.format(max_product))
    ax1.set_ylabel('H')
    ax1.grid(True)
    ax1.legend(loc='upper right')
    legends = ax1.get_legend()
    legends.legendHandles[0].set_color('xkcd:azure')
    legends.legendHandles[1].set_color('xkcd:orangered')
    ax2 = ax1.twinx()
    ax2.plot(index, count_min_len, 'b-', label='{}bp'.format(min_product))
    ax2.plot(index, count_max_len, 'r-', label='{}bp'.format(max_product))
    ax2.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0))
    ax2.set_ylabel('Resolution(% of {})'.format(rows))
    ax2.grid(True)
    ax2.legend(loc='upper left')
    plt.savefig(out+'.pdf')
    plt.savefig(out+'.png')
    # plt.show()
    with open(out+'-Resolution.tsv', 'w') as _:
        _.write('Base\tResolution(window={})\n'.format(min_product))
        for base, resolution in enumerate(count_min_len):
            _.write('{}\t{:.2f}\n'.format(base, resolution))

    return (count_min_len, count_max_len, shannon_index1, shannon_index2,
            max_shannon_index, index)


#@profile
def validate(query_file, db_file, n_seqs, min_len, min_covrage,
             max_mismatch):
    """
    Do BLAST. Parse BLAST result. Return Dict[str, Dict[str, Any]].
    """
    # build blast db
    run('makeblastdb -in {} -dbtype nucl'.format(db_file), shell=True,
        stdout=open('makeblastdb.log', 'w'))
    # blast
    blast_result_file = 'BlastResult.xml'
    cmd = nb(num_threads=cpu_count(),
             query=query_file,
             db=db_file,
             task='blastn',
             evalue=1e-5,
             max_hsps=1,
             max_target_seqs=n_seqs,
             outfmt=5,
             out=blast_result_file)
    stdout, stderr = cmd()
    # minium match bases * score per base(2)
    min_bitscore_raw = (min_len - max_mismatch)*2
    blast_result: Dict[int, Dict[str, Any]] = dict()
    for query in SearchIO.parse(blast_result_file, 'blast-xml'):
        index = int(query.id.split('-')[0])
        if len(query) == 0:
            blast_result[index] = {'coverage': 0, 'sum_bitscore': 0,
                                   'avg_mid_loc': 0}
            continue
        sum_bitscore_raw = 0
        good_hits = 0
        start = 0
        for hit in query:
            hsp_bitscore_raw = hit[0].bitscore_raw
            if hsp_bitscore_raw >= min_bitscore_raw:
                sum_bitscore_raw += hsp_bitscore_raw
                good_hits += 1
                start += sum(hit[0].hit_range) / 2
        coverage = good_hits/n_seqs
        if coverage >= min_covrage:
            # get sequence unique index for choose
            # Notice that here it use bitscore_raw instead of bitscore
            blast_result[index] = {'coverage': coverage,
                                   'sum_bitscore': sum_bitscore_raw,
                                   'avg_mid_loc': start/n_seqs}
    return blast_result


#@profile
def parse_args():
    arg = argparse.ArgumentParser(description=main.__doc__)
    arg.add_argument('input', help='input alignment file')
    arg.add_argument('-a', '--ambiguous_base_n', type=int, default=2,
                     help='number of ambiguous bases')
    arg.add_argument('-c', '--cutoff', type=float, default=0.8,
                     help='minium percent to keep base')
    arg.add_argument('-g', '--gap_cutoff', type=float, default=0.5,
                     help='maximum percent for gap to cutoff')
    arg.add_argument('-pmin', '--min_primer', type=int, default=24,
                     help='minimum primer length')
    arg.add_argument('-pmax', '--max_primer', type=int, default=25,
                     help='maximum primer length')
    arg.add_argument('-m', '--mismatch', type=int, default=2,
                     help='maximum mismatch bases in primer')
    arg.add_argument('-o', '--out', help='output name prefix')
    arg.add_argument('-r', '--resolution', type=float, default=0.6,
                     help='minium resolution')
    arg.add_argument('-tmin', '--min_product', type=int, default=380,
                     help='minimum product length(include primer)')
    arg.add_argument('-tmax', '--max_product', type=int, default=480,
                     help='maximum product length(include primer)')
    arg.add_argument('-w', '--window', type=int, default=1,
                     help='window size')
    # arg.print_help()
    return arg.parse_args()


#@profile
def main():
    """
    Automatic design primer for DNA barcode.
    """
    start = timer()
    arg = parse_args()
    if arg.out is None:
        arg.out = os.path.basename(arg.input)
        arg.out = arg.out.split('.')[0]

    # read from fasta, generate new fasta for makeblastdb
    name, alignment, db_file = prepare(arg.input)
    rows, columns = alignment.shape


    # generate consensus
    base_cumulative_frequency = count_base(alignment, rows, columns)
    consensus = generate_consensus(base_cumulative_frequency, arg.cutoff,
                                   arg.gap_cutoff, rows, columns,
                                   arg.out+'.consensus.fastq')

    # count resolution
    (seq_count_min_len, seq_count_max_len,
     H1, H2, max_H, index) = count_and_draw(alignment, consensus, arg)

    # exit if resolution lower than given threshold.
    assert max(seq_count_max_len) > arg.resolution, (
        """
The highest resolution of given fragment is {:.2f}, which is lower than
given resolution threshold({:.2f}). Please try to use longer fragment or
lower resolution options.
""".format(max(seq_count_max_len), arg.resolution))

    good_region = [False] * columns
    n = arg.max_product - arg.min_product + arg.max_primer
    # lower bound, min_prodcut with max_primer
    # upper bound, max_prodcut with min_primer
    n2 = arg.max_product + arg.min_primer
    for i, j, k in zip(index, seq_count_min_len, seq_count_max_len):
        if j >= arg.resolution:
            for _ in range(i-n, i):
                good_region[_] = True
            for _ in range(i+arg.min_product, i+arg.min_product+arg.max_primer):
                good_region[_] = True
        elif k >= arg.resolution:
            for _ in range(i-arg.min_primer, i):
                good_region[_] = True
            for _ in range(i+arg.max_product, i+n2):
                good_region[_] = True
                # +1 or not???????
    # find candidate
    continuous = find_continuous(consensus, good_region, arg.min_primer)
    primer_candidate = find_primer(continuous, rows, arg.min_primer,
                                   arg.max_primer, arg.ambiguous_base_n)
    assert len(primer_candidate) != 0, (
        'Primer not found! Try to loose options.')
    candidate_file = arg.out + '.candidate.fasta'
    with open(candidate_file, 'w') as _:
        for item in primer_candidate:
            item.write(_, 'fasta')

    # validate
    result = validate(candidate_file, db_file, rows, arg.min_primer,
                      arg.cutoff, arg.mismatch)
    primer_file = '{}-{}_covrage-{}bp_mismatch.fastq'.format(
        arg.out, arg.cutoff, arg.mismatch)
    count = 0
    with open(primer_file, 'w') as out:
        for primer in primer_candidate:
            i = primer.index
            if i in result:
                count += 1
                primer.coverage = result[i]['coverage']
                primer.sum_bitscore = result[i]['sum_bitscore']
                primer.avg_mid_loc = result[i]['avg_mid_loc']
                primer.write(out, 'fastq')

    print('Found {} primers.'.format(count))
    print('Primer ID format:')
    print('name-Tm-Samples-BLAST_Coverage-Bitscore-AvgMidLocation')
    end = timer()
    print('Cost {:.3f} seconds.'.format(end-start))


if __name__ == '__main__':
    main()
