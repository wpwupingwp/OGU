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

from matplotlib import pyplot as plt
from matplotlib import ticker as mtick
import matplotlib
matplotlib.rcParams['lines.linewidth'] = 1.5
matplotlib.rcParams['axes.linewidth'] = 1.5
matplotlib.rcParams['axes.labelsize'] = 16
matplotlib.rcParams['axes.titlesize'] = 25
matplotlib.rcParams['font.size'] = 16
matplotlib.rcParams['axes.facecolor'] = '#888888'
matplotlib.rcParams['figure.figsize'] = 16, 9


class PrimerInfo:
    def __init__(self, index=None, start=None, end=None, tm=None,
                 coverage=None, bitscore=None, avg_mid_location=None,
                 detail=None):
        self.id = index
        self.start = start
        self.end = end
        self.tm = tm
        self.coverage = coverage
        selft.bitscore = bitscore
        self.avg_mid_location = avg_mid_location
        self.detail = detail

    def __str__(self):
        return ('No.{}-Start:{}-End:{}-Tm:{}℃-Coverage:{:.2%}-Bitscore'
                ':{}-AvgMidLocation:{}')


def prepare(fasta):
    """
    Given fasta format alignment filename, return a numpy array for sequence:
    List[ID, Sequence].
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
    if len(set(length_check)) != 1:
        raise ValueError(
            'Alignment does not have uniform width, please check again !')

    # Convert List to numpy array.
    # order 'F' is a bit faster than 'C'
    # new = np.hstack((name, seq)) -> is slower
    name = np.array([[i[0]] for i in data], dtype=np.bytes_)
    sequence = np.array([list(i[1]) for i in data], dtype=np.bytes_, order='F')
    return name, sequence, no_gap


def count(alignment, rows, columns):
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


def get_quality_string(data: List[float], rows: int):
    # https://en.wikipedia.org/wiki/FASTQ_format
    quality = ('''!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ'''
               '''[\]^_`abcdefghijklmnopqrstuvwxyz{|}~''')
    quality_dict = {i: j for i, j in enumerate(quality)}
    max_q = len(quality)
    factor = max_q/rows
    # use min to avoid KeyError
    quality_value = [min(max_q, int(i*factor))-1 for i in data]
    quality_string = [quality_dict[i] for i in quality_value]
    return  ''.join(quality_string)


def generate_consensus(base_cumulative_frequency, cutoff, gap_cutoff,
                       rows, columns):
    """
    Given base count info, return List[index, base, quality]
    and List[List[str, str, str]] for writing conesensus.
    """
    # directly use np.unique result
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
    consensus = [0, ]
    consensus.append(''.join([i[1] for i in most]))
    quality = [i[2] for i in most]
    consensus.append(get_quality_string(quality))
    # try to use class
    consensus.append(PrimerInfo(start=1, end=len(most), coverage=cutoff,
                                avg_mid_location=len(most)/2)
    # write require List[List[]]
    return most, [consensus, ]


# to be continued
def find_continuous(consensus, min_len):
    """
    Given List[location, base, count]
    Return continuous fragments list:
    List[List[location, base, count]]
    """
    continuous: List[List[int, str, float]] = list()
    fragment = list()
    skip = ('N', '-', '*')
    for base in consensus:
        if base[1] in skip:
            if len(fragment) >= min_len:
                continuous.append(fragment)
                fragment = list()
            continue
        fragment.append(base)
    return continuous


def find_primer(continuous, most, min_len, max_len, ambiguous_base_n):
    """
    Find suitable primer in given List[List[int, str, float]]
    return List[List[int, str, List[float], Dict[str, Any]]]
    """
    poly = re.compile(r'([ATCG])\1\1\1\1')
    ambiguous_base = re.compile(r'[^ATCG]')
    tandem = re.compile(r'([ATCG]{2})\1\1\1\1')

    def is_good_primer(primer):
        # ref1. http://www.premierbiosoft.com/tech_notes/PCR_Primer_Design.html
        seq = ''.join([i[1] for i in primer])
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

    primer: List[List[int, str, Dict[str, Any]]] = list()
    continuous = [i for i in continuous if len(i) >= min_len]
    n = 1
    for fragment in continuous:
        len_fragment = len(fragment)
        for begin in range(len_fragment-max_len):
            for p_len in range(min_len, max_len):
                seq = fragment[begin:(begin+p_len)]
                good_primer, tm, detail = is_good_primer(seq)
                if good_primer:
                    start = seq[0][0]
                    end = seq[-1][0]
                    sequence = ''.join([i[1] for i in seq])
                    quality = [i[2] for i in seq]
                    primer.append([n, sequence, get_quality_string(quality),
                                   PrimerInfo(start=start, end=end, tm=tm,])
                    n += 1
                else:
                    continue
    return primer


def unique_sequence_count(data, window):
    rows, columns = data.shape
    # Different count
    C = list()
    factor = 100/rows
    for i in range(columns-window):
        cut = data[:, i:(i+window)]
        # uniqe array, count line*times
        _, count = np.unique(cut, return_counts=True, axis=0)
        C.append(len(count)*factor)
    return C


def shannon_diversity_index(data, rows, columns, sequence_count_result, window,
                            only_atcg=True, with_n=False, with_gap=False,
                            out='out.png'):
    """http://www.tiem.utk.edu/~gross/bioed/bealsmodules/shannonDI.html
    """
    # only_atcg: only consider ATCG 4 kinds of bases
    # with_n: consider N as the fifth kind of base
    # with_gap: consider N as the fifth kind of base and gap as the sixth kind
    # of base

    # split shape and data
    data = np.array(data)
    if with_gap:
        new_data = data[:, 0:6]
    elif with_n:
        new_data = data[:, 0:5]
    elif only_atcg:
        new_data = data[:, 0:4]
    # Shannon Index
    H = list()
    # Sum_all/max_h
    size = list()
    # max_h = -1*((1/len(new_data[0]))*log2(1/(len(new_data[0]))))*len(
    #     new_data[0])
    # dot size
    max_size = 50
    for column in new_data:
        # sum_all equals sum of letters considered rather than original rows
        sum_all = sum(column)
        size.append(sum_all/rows*max_size)
        h = 0
        for i in column:
            if i == 0:
                continue
            p_i = i / sum_all
            log2_p_i = log2(p_i)
            h += log2_p_i*p_i
        H.append(-1*h)
    # plt.style.use('ggplot')
    fig, ax1 = plt.subplots()
    plt.title('Shannon Diversity Index & Resolution(window={})'.format(window))
    plt.xlabel('Base')
    plt.xticks(range(0, columns, int(columns/10)))
    # ax1.plot((0, columns), (max_h, max_h), 'r--', label='Max H')
    # c=List for different color, s=size for different size
    ax1.scatter(range(columns), H, c=H, cmap='GnBu', s=size)
    ax1.set_ylabel('H')
    ax1.grid(True)
    ax2 = ax1.twinx()
    ax2.plot(sequence_count_result, 'r-', alpha=0.8)
    ax2.yaxis.set_major_formatter(mtick.PercentFormatter())
    ax2.set_ylabel('Resolution(% of {})'.format(rows))
    ax2.grid(True)
    # plt.legend(loc=1, frameon=False)
    plt.savefig(out+'.pdf')
    plt.savefig(out+'.png')
    # plt.show()
    with open(out+'-Resolution.tsv', 'w') as _:
        _.write('Base\tResolution(window={})\n'.format(window))
        for base, resolution in enumerate(sequence_count_result):
            _.write('{}\t{:.2f}\n'.format(base, resolution))


def validate(query_file, db_file, n_seqs, min_len, min_covrage,
             max_mismatch):
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
    blast_result = [['ID', 'Hits', 'Sum_Bitscore_raw'], ]
    blast_result.append(['All', n_seqs, min_len])
    for query in SearchIO.parse(blast_result_file, 'blast-xml'):
        if len(query) == 0:
            blast_result.append([query.id, 0, 0])
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
        blast_result.append([query.id, good_hits/n_seqs, sum_bitscore_raw,
                             start/n_seqs])
    # validate
    # validate_result = [['ID', 'Hits', 'Sum_Bitscore_raw', 'Seq'], ]
    validate_result = list()
    for record in blast_result[2:]:
        if record[1] >= min_covrage:
            validate_result.append(record)
    validate_result.sort(key=lambda x: x[1], reverse=True)
    return validate_result


def write_to_file(data, rows, output, name, file_format):
    """
    Given List[List[int, str, str, Dict[str, Any]]
    Write fasta or fastq format file.
    """
    out = open(output, 'w')
    for item in data:
        for key, value in data[-1].items():
        sequence_id = 'No.{}-{}'.format(item[0], 
        if file_format == 'fastq':
            out.write('@{}-{}-{}-{:.3f}-{}\n'.format(name, start, end, tm, rows))
            out.write(seq+'\n')
            out.write('+\n')
            out.write(''.join(qual_character)+'\n')
        elif file_format == 'fasta':
            out.write('>{}-{}-{}-{:.3f}-{}\n'.format(name, start, end, tm, rows))
            out.write(seq+'\n')
    return output


def parse_args():
    arg = argparse.ArgumentParser(description=main.__doc__)
    arg.add_argument('input', help='input alignment file')
    arg.add_argument('-a', '--ambiguous_base_n', type=int, default=2,
                     help='number of ambiguous bases')
    arg.add_argument('-c', '--cutoff', type=float, default=0.9,
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
    # todo
    arg.add_argument('-r', '--resolution', help='minium resolution')
    arg.add_argument('-tmin', '--min_template', type=int, default=350,
                     help='minimum template length')
    arg.add_argument('-tmax', '--max_template', type=int, default=450,
                     help='maximum template length')
    # arg.print_help()
    return arg.parse_args()


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
    print(rows, columns)

    # generate consensus
    base_cumulative_frequency = count(alignment, rows, columns)
    consensus, consensus_for_write = generate_consensus(
        base_cumulative_frequency, arg.cutoff, arg.gap_cutoff, rows, columns)

    # find candidate
    continuous = find_continuous(consensus, arg.min_primer)
    primer_candidate = find_primer(continuous, consensus, arg.min_primer,
                                   arg.max_primer, arg.ambiguous_base_n)
    if len(primer_candidate) == 0:
        raise ValueError('Primer not found! Try to loose restriction.')
    candidate_file = write_to_file(primer_candidate, rows,
                                   arg.out+'.candidate.fasta', arg.out,
                                   'fasta')
    candidate_file_fastq = write_to_file(primer_candidate, rows,
                                         arg.out+'.candidate.fastq', arg.out,
                                         'fastq')

    # validate
    primer_info = validate(candidate_file, db_file, rows, arg.min_primer,
                           arg.cutoff, arg.mismatch)
    primer_info_dict = {i[0]: i[1:] for i in primer_info}
    primer_file = '{}-{}_covrage-{}bp_mismatch.fastq'.format(
        arg.out, arg.cutoff, arg.mismatch)
    # write
    with open(primer_file, 'w') as out:
        for seq in SeqIO.parse(candidate_file_fastq, 'fastq'):
            if seq.id in primer_info_dict:
                short_id = seq.id.split('-')
                short_id = '-'.join([short_id[0], short_id[-2], short_id[-1]])
                # name-Tm-Samples-BLAST_Coverage-Bitscore-AvgMidLocation
                seq.id = '{}-{:.2%}-{}-{:.2f}'.format(
                    short_id, *primer_info_dict[seq.id])
                seq.description = ''
                SeqIO.write(seq, out, 'fastq')
    sequence_count_result = unique_sequence_count(alignment,
                                                  window=arg.min_template)
    shannon_diversity_index(base_cumulative_frequency, rows, columns,
                            sequence_count_result, window=arg.min_template,
                            only_atcg=True, out=arg.out)

    # write consensus
    write_to_file(consensus_for_write, rows, arg.out+'.consensus.fastq',
                  arg.out, 'fastq')
    print('Found {} primers.'.format(len(primer_info)))
    print('Primer ID format:')
    print('name-Tm-Samples-BLAST_Coverage-Bitscore-AvgMidLocation')
    end = timer()
    print('Cost {:.3f} seconds.'.format(end-start))


if __name__ == '__main__':
    main()
