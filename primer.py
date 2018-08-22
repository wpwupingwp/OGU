#!/usr/bin/python3

import argparse
import json
import numpy as np
import os
import primer3
import re
from Bio import Phylo, SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline as Blast
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import defaultdict
from glob import glob
from matplotlib import pyplot as plt
from matplotlib import rcParams
from subprocess import run
from tempfile import NamedTemporaryFile as Tmp
from timeit import default_timer as timer

rcParams['lines.linewidth'] = 1.5
rcParams['axes.linewidth'] = 1.5
rcParams['axes.labelsize'] = 16
rcParams['axes.titlesize'] = 25
rcParams['font.size'] = 16


class PrimerWithInfo(SeqRecord):
    def __init__(self, seq='', quality='', start=0, coverage=0, avg_bitscore=0,
                 mid_loc=None, avg_mismatch=0, detail=0,
                 reverse_complement=False):
        # store str
        super().__init__(Seq(seq.upper()))
        self.sequence = str(self.seq)

        """
        primer3.setGlobals seems have no effect on calcTm, so I have to replace
        all ambiguous base to G to get an approximate value. Othervise calcTm()
        will generate -99999 if there is ambiguous base.
        """
        self.g_seq = re.sub(r'[^ATCG]', 'G', self.sequence)
        self.quality = self.letter_annotations['solexa_quality'] = quality
        self.start = self.annotations['start'] = start
        self.tm = self.annotations['tm'] = primer3.calcTm(self.g_seq)
        self.coverage = self.annotations['coverage'] = coverage
        self.avg_bitscore = self.annotations['avg_bitscore'] = avg_bitscore
        self.mid_loc = self.annotations['mid_loc'] = mid_loc
        self.avg_mid_loc = 0
        self.avg_mismatch = self.annotations['avg_mismatch'] = avg_mismatch
        self.detail = self.annotations['detail'] = detail
        self.end = self.annotations['end'] = start + self.__len__() - 1
        self.is_reverse_complement = reverse_complement
        self.description = ''
        self.hairpin_tm = 0
        self.homodimer_tm = 0
        self.update_id()

    def __getitem__(self, i):
        if isinstance(i, int):
            i = slice(i, i+1)
        if isinstance(i, slice):
            if self.seq is None:
                raise ValueError('Empty sequence')
            answer = PrimerWithInfo(seq=str(self.seq[i]),
                                    quality=self.quality[i])
            answer.annotations = dict(self.annotations.items())
            return answer

    def is_good_primer(self):
        poly = re.compile(r'([ATCG])\1\1\1\1')
        tandem = re.compile(r'([ATCG]{2})\1\1\1\1')
        # ref1. http://www.premierbiosoft.com/tech_notes/PCR_Primer_Design.html
        if re.search(poly, self.g_seq) is not None:
            self.detail = 'Poly(NNNNN) structure found'
            return False
        if re.search(tandem, self.g_seq) is not None:
            self.detail = 'Tandom(NN*5) exist'
            return False
        self.hairpin_tm = primer3.calcHairpinTm(self.g_seq)
        self.homodimer_tm = primer3.calcHomodimerTm(self.g_seq)
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
        table = str.maketrans('ATCGRMYKHBDV', 'TAGCYKRMDVHB')
        new_seq = str.translate(self.sequence, table)
        new_seq = new_seq[::-1]
        new_quality = self.quality[::-1]
        # try to simplify??
        return PrimerWithInfo(seq=new_seq, quality=new_quality,
                              start=self.start, coverage=self.coverage,
                              avg_bitscore=self.avg_bitscore,
                              mid_loc=self.mid_loc,
                              detail=self.detail)

    def update_id(self):
        self.end = self.annotations['end'] = self.start + self.__len__() - 1
        if self.mid_loc is not None and len(self.mid_loc) != 0:
            self.avg_mid_loc = np.average(list(self.mid_loc.values()))
        self.id = ('AvgMidLocation({:.0f})-Tm({:.2f})-Coverage({:.2%})-'
                   'AvgBitScore({:.2f})-Start({})-End({})'.format(
                       self.avg_mid_loc, self.tm, self.coverage,
                       self.avg_bitscore, self.start, self.end))


class Pair:
    # save memory
    __slots__ = ['left', 'right', 'delta_tm', 'coverage', 'start', 'end',
                 'resolution', 'tree_value', 'entropy', 'have_heterodimer',
                 'heterodimer_tm', 'pi', 'score', 'length']

    def __init__(self, left, right, alignment):
        rows, columns = alignment.shape
        self.left = left
        if not right.is_reverse_complement:
            self.right = right.reverse_complement()
        else:
            self.right = right
        self.delta_tm = abs(self.left.tm - self.right.tm)
        a = len(self.left)/2
        b = len(self.right)/2
        lengths = dict()
        common = self.left.mid_loc.keys() & self.right.mid_loc.keys()
        for key in common:
            length = (self.right.mid_loc[key]-b) - (self.left.mid_loc[key]+a)
            # omit negative length
            if length > 0:
                lengths[key] = int(length)
        self.length = lengths
        self.left.mid_loc = {i: j for i, j in self.left.mid_loc.items() if
                             (i in lengths)}
        self.right.mid_loc = {i: j for i, j in self.right.mid_loc.items() if
                              (i in lengths)}
        self.left.coverage = len(self.left.mid_loc) / rows
        # recalculate coverage due to dropping some records
        self.right.coverage = len(self.right.mid_loc) / rows
        self.left.update_id()
        self.right.update_id()
        self.coverage = min(self.left.coverage, self.right.coverage)
        # pairs use mid_loc from BLAST as start/end
        self.start = int(self.left.avg_mid_loc)
        self.end = int(self.right.avg_mid_loc)
        self.resolution = 0
        self.tree_value = 0.0
        self.entropy = 0.0
        self.have_heterodimer = False
        self.resolution = 0
        self.entropy = 0
        self.pi = 0.0
        self.heterodimer_tm = primer3.calcHeterodimer(self.left.g_seq,
                                                      self.right.g_seq).tm
        if max(self.heterodimer_tm, self.left.tm,
               self.right.tm) == self.heterodimer_tm:
            self.have_heterodimer = True
        else:
            self.have_heterodimer = False
        self.get_score()

    def __str__(self):
        return (
            'Pair(score={:.2f}, product={:.0f}, start={}, end={}, left={}, '
            'right={}, resolution={:.2%}, coverage={:.2%}, delta_tm={:.2f}, '
            'have_heterodimer={})'.format(
                self.score, np.average(list(self.length.values())), self.start,
                self.end, self.left.seq, self.right.seq, self.resolution,
                self.coverage, self.delta_tm, self.have_heterodimer))

    def get_score(self):
        self.score = (np.average(list(self.length.values()))*0.5
                      + self.coverage*200
                      + len(self.left)*10
                      + len(self.right)*10
                      + self.resolution*100
                      + self.tree_value*100 + self.entropy*5
                      - int(self.have_heterodimer)*10
                      - self.delta_tm*5 - self.left.avg_mismatch*10
                      - self.right.avg_mismatch*10)

    def add_info(self, alignment):
        # include end base, use alignment loc for slice
        (self.resolution, self.entropy, self.pi,
         self.tree_value) = get_resolution(alignment, self.left.start,
                                           self.right.end+1)
        self.get_score()
        return self


class BlastResult():
    __slots = ['query_id', 'hit_id', 'query_seq', 'ident_num', 'mismatch_num',
               'bitscore_raw', 'query_start', 'query_end', 'hit_start',
               'hit_end']

    def __init__(self, line):
        record = line.strip().split('\t')
        self.query_id, self.hit_id, self.query_seq = record[0:3]
        (self.ident_num, self.mismatch_num, self.bitscore_raw,
         self.query_start, self.query_end, self.hit_start,
         self.hit_end) = [int(i) for i in record[3:]]


def parse_args():
    arg = argparse.ArgumentParser(
        description=main.__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    arg.add_argument('input', help='input alignment file')
    options = arg.add_argument_group()
    options.add_argument('-a', '--ambiguous_base_n', type=int, default=4,
                         help='number of ambiguous bases')
    options.add_argument('-c', '--coverage', type=float, default=0.6,
                         help='minium coverage of base and primer')
    options.add_argument('-f', '--fast', action='store_true', default=False,
                         help='faster evaluate variance by omit tree_value')
    options.add_argument('-j', '--json', help='configuration json file')
    options.add_argument('-m', '--mismatch', type=int, default=4,
                         help='maximum mismatch bases in primer')
    options.add_argument('-pmin', '--min_primer', type=int, default=18,
                         help='minimum primer length')
    options.add_argument('-pmax', '--max_primer', type=int, default=24,
                         help='maximum primer length')
    options.add_argument('-r', '--resolution', type=float, default=0.5,
                         help='minium resolution')
    options.add_argument('-s', '--step', type=int, default=50,
                         help='step size')
    options.add_argument('-t', '--top_n', type=int, default=1,
                         help='keep n primers for each high varient region')
    options.add_argument('-tmin', '--min_product', type=int, default=300,
                         help='minimum product length(include primer)')
    options.add_argument('-tmax', '--max_product', type=int, default=500,
                         help='maximum product length(include primer)')
    # arg.print_help()
    parsed = arg.parse_args()
    parsed.out = os.path.basename(parsed.input)
    parsed.out = parsed.out.split('.')[:-1]
    parsed.out = '.'.join(parsed.out)
    # overwrite options by given json
    if parsed.json is not None:
        with open(parsed.json, 'r') as _:
            config = json.load(_)
        n_arg = argparse.Namespace(**config)
        n_arg.input = parsed.input
        n_arg.out = parsed.out
        return n_arg
    else:
        return parsed


def prepare(fasta):
    """
    Given fasta format alignment filename, return a numpy array for sequence:
    Generate fasta file without gap for makeblastdb, return file name.
    """
    no_gap = Tmp('wt')
    data = list()
    record = ['id', 'sequence']
    with open(fasta, 'r') as raw:
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
    organize_no_gap = Tmp('wt', delete=False)
    # try to avoid makeblastdb error
    SeqIO.convert(no_gap.name, 'fasta', organize_no_gap.name, 'fasta')
    no_gap.close()
    return name, sequence, organize_no_gap.name


def count_base(alignment, rows, columns):
    """
    Given alignment numpy array, count cumulative frequency of base in each
    column (consider ambiguous base and "N", "-" and "?", otherwise omit).
    Return List[List[float, float, float, float, float, float, float]] for
    [A, T, C, G, N, GAP, OTHER].
    """
    frequency = list()
    for index in range(columns):
        base, counts = np.unique(alignment[:, [index]], return_counts=True)
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


def get_quality(data, rows):
    # use fastq-illumina format
    max_q = 62
    factor = max_q/rows
    # use min to avoid KeyError
    quality_value = [min(max_q, int(i*factor))-1 for i in data]
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
        d_ij = np.sum(subalign[i] != subalign[(i+1):])
        sum_d_ij += d_ij
    pi = (2 / (n*(n-1)) * sum_d_ij) / m
    # tree value
    aln_file = '{}-{}.aln.tmp'.format(start, end)

    def clean():
        for i in glob(aln_file+'*'):
            os.remove(i)
    if not fast:
        with open(aln_file, 'wb') as aln:
            for index, row in enumerate(alignment[:, start:end]):
                aln.write(b'>'+str(index).encode('utf-8')+b'\n'+b''.join(
                    row)+b'\n')
        iqtree = run('iqtree -s {} -m JC -fast -czb'.format(aln_file),
                     stdout=Tmp('wt'), shell=True)
        # just return 0 if there is error
        if iqtree.returncode != 0:
            print('Cannot get tree_value of {}-{}!'.format(start, end))
            clean()
            return resolution, entropy, pi, 0
        tree = Phylo.read(aln.name+'.treefile', 'newick')
        # skip the first empty node
        internals = tree.get_nonterminals()[1:]
        clean()
        tree_value = len(internals) / len(tree.get_terminals())
    else:
        tree_value = 0
    return resolution, entropy, pi, tree_value


def generate_consensus(base_cumulative_frequency, coverage_percent,
                       rows, columns, output):
    """
    Given base count info, return List[index, base, quality]
    and List[List[str, str, str, PrimerInfo]] for writing conesensus.
    return PrimerWithInfo
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
    most = list()
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
    quality = [i[2] for i in most]
    consensus = PrimerWithInfo(start=1, seq=''.join([i[1] for i in most]),
                               quality=get_quality(quality, rows))
    SeqIO.write(consensus, output, 'fastq')
    return consensus


def get_good_region(index, seq_count, arg):
    # return loose region, final product may violate product length
    # restriction
    n = arg.max_product - arg.min_product
    good_region = set()
    for i, j in zip(index, seq_count):
        if j >= arg.resolution:
            good_region.update(range(i-arg.max_primer, i-n))
            good_region.update(range(i+arg.min_product,
                                     i-arg.max_primer+arg.max_product))
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
            if (index-start) >= min_len:
                consensus.features.append(SeqFeature(FeatureLocation(
                    start, index), type='continuous', strand=1))
            start = index + 1
    return consensus


def find_primer(consensus, min_len, max_len):
    """
    Find suitable primer in given consensus with features labeled as candidate
    primer, return List[PrimerWithInfo], consensus
    """
    primers = list()
    # skip good_region
    continuous = consensus.features
    for feature in continuous:
        fragment = feature.extract(consensus)
        len_fragment = len(fragment)
        for begin in range(len_fragment-max_len):
            for p_len in range(min_len, max_len+1):
                start = feature.location.start + begin
                primer = consensus[start:start+p_len]
                if primer.is_good_primer():
                    consensus.features.append(SeqFeature(
                        FeatureLocation(start, start+p_len),
                        type='primer', strand=1))
                    primer.start = start
                    primer.update_id()
                    primers.append(primer)
    return primers, consensus


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
    max_product = arg.max_product
    step = arg.step
    out = arg.out
    # R, H, Pi, T : count, normalized entropy, Pi, tree value
    R = list()
    H = list()
    Pi = list()
    T = list()
    max_H = np.log2(rows)
    index = list()
    max_plus = max_product - min_primer * 2
    max_range = columns - max_product
    for i in range(0, max_range, step):
        # skip gap
        # if consensus.sequence[i] in ('-', 'N'):
        #     continue
        # exclude primer sequence
        resolution, entropy, pi, tree_value = get_resolution(
            alignment, i, i+max_plus, arg.fast)
        R.append(resolution)
        H.append(entropy/max_H)
        Pi.append(pi)
        T.append(tree_value)
        index.append(i)

    plt.style.use('seaborn-colorblind')
    fig, ax1 = plt.subplots(figsize=(15+len(index)//5000, 10))
    plt.title('Resolution(window={} bp, step={} bp)\n'.format(
        max_product, step))
    plt.xlabel('Base')
    # plt.xticks(np.linspace(0, max_range, 21))
    if not arg.fast:
        ax1.set_ylabel('Normalized Shannon Index / Resolution / TreeValue')
        ax1.plot(index, T, label='TreeValue', alpha=0.8)
    else:
        ax1.set_ylabel('Normalized Shannon Index / Resolution')

    ax1.plot(index, H, label='Shannon Index', alpha=0.8)
    ax1.plot(index, R, label='Resolution', alpha=0.8)
    ax1.legend(loc='lower left')
    ax1.yaxis.set_ticks(np.linspace(0, 1, num=11))
    ax2 = ax1.twinx()
    ax2.plot(index, Pi, 'k-', label=r'$\pi$', alpha=0.8)
    ax2.set_ylabel(r'$\pi$', rotation=-90, labelpad=20)
    _ = round(np.log10(max(Pi)))
    ax2.yaxis.set_ticks(np.linspace(0, 10**_, num=11))
    ax2.legend(loc='upper right')
    # plt.yscale('log')
    plt.savefig(out+'.pdf')
    plt.savefig(out+'.png')
    # plt.show()
    with open(out+'-Resolution.tsv', 'w') as _:
        _.write('Index,R,H,Pi,T\n')
        for i, r, h, pi, t in zip(index, R, H, Pi, T):
            _.write('{},{:.2f},{:.2f},{:.2f},{:.2f}\n'.format(i, r, h, pi, t))
    return R, H, Pi, T, index


def parse_blast_tab(filename):
    query = list()
    with open(filename, 'r') as raw:
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
    query_file = arg.out + '.candidate.fasta'
    # SeqIO.write fasta file directly is prohibited. have to write fastq at
    with open(query_file+'.fastq', 'w') as _:
        SeqIO.write(primer_candidate, _, 'fastq')
    SeqIO.convert(query_file+'.fastq', 'fastq', query_file, 'fasta')
    # build blast db
    run('makeblastdb -in {} -dbtype nucl'.format(db_file), shell=True,
        stdout=Tmp('wt'))
    # blast
    blast_result_file = Tmp('wt')
    fmt = 'qseqid sseqid qseq nident mismatch score qstart qend sstart send'
    cmd = Blast(num_threads=len(os.sched_getaffinity(0)),
                query=query_file,
                db=db_file,
                task='blastn-short',
                evalue=1e-2,
                max_hsps=1,
                outfmt='"7 {}"'.format(fmt),
                out=blast_result_file.name)
    stdout, stderr = cmd()
    blast_result = dict()
    # because SearchIO.parse is slow, use parse_blast_result()
    for query in parse_blast_tab(blast_result_file.name):
        if len(query) == 0:
            continue
        sum_bitscore_raw = 0
        sum_mismatch = 0
        good_hits = 0
        mid_loc = dict()
        for hit in query:
            min_positive = len(hit.query_seq) - arg.mismatch
            hsp_bitscore_raw = hit.bitscore_raw
            positive = hit.ident_num
            mismatch = hit.mismatch_num
            loc = np.average([hit.hit_start, hit.hit_end])
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
                'avg_bitscore': sum_bitscore_raw/good_hits,
                'avg_mismatch': sum_mismatch/good_hits,
                'mid_loc': mid_loc}
    primer_verified = list()
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
    blast_result_file.close()
    # clean makeblastdb files
    for i in glob(db_file+'*'):
        os.remove(i)
    return primer_verified


def pick_pair(primers, alignment, arg):
    pairs = list()
    for n_left, left in enumerate(primers):
        # convert mid_loc to 5' location
        location = left.avg_mid_loc - len(left) / 2
        begin = location + arg.min_product
        # fragment plus one primer = max_product length
        end = location + arg.max_product - len(left)
        cluster = list()
        for right in primers[(n_left+1):]:
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
    assert len(pairs) != 0, 'Primers pairs not found!'
    # remove close located primers
    less_pairs = list()
    cluster = [pairs[0], ]
    for index in range(1, len(pairs)):
        if pairs[index].start - pairs[index-1].start < arg.min_primer:
            cluster.append(pairs[index])
        else:
            cluster.sort(key=lambda x: x.score, reverse=True)
            less_pairs.extend(cluster[:arg.top_n])
            cluster.clear()
    cluster.sort(key=lambda x: x.score, reverse=True)
    less_pairs.extend(cluster[:arg.top_n])
    good_pairs = list()
    for i in less_pairs:
        i.add_info(alignment)
        if i.resolution >= arg.resolution:
            good_pairs.append(i)
    good_pairs.sort(key=lambda x: x.score, reverse=True)
    return good_pairs


def main():
    """
    Automatic design primer for DNA barcode.
    """
    summary = 'Summary.csv'
    start = timer()
    arg = parse_args()
    print('Write configuration into json')
    with open(arg.input+'.json', 'w') as out:
        json.dump(vars(arg), out, indent=4, sort_keys=True)
    # read from fasta, generate new fasta for makeblastdb
    name, alignment, db_file = prepare(arg.input)
    rows, columns = alignment.shape
    assert rows >= 4, 'Too few sequence in alignment (less than 4)!'
    # generate consensus
    base_cumulative_frequency = count_base(alignment, rows, columns)
    consensus = generate_consensus(base_cumulative_frequency, arg.coverage,
                                   rows, columns, arg.out+'.consensus.fastq')
    max_count, max_H, max_Pi, max_T = get_resolution(alignment, 0, columns)
    n_gap = sum([i[5] for i in base_cumulative_frequency])
    gap_ratio = n_gap / rows / columns
    if not os.path.exists(summary):
        with open(summary, 'w') as s:
            s.write('Name,Sequences,Length,GapRatio,ObservedResolution,'
                    'TreeValue,ShannonIndex,Pi\n')
            s.write('{},{},{},{:.2%},{:.6f},{:.6f},{:.6f},{:.6f}\n'.format(
                os.path.basename(arg.input), rows, columns, gap_ratio,
                max_count, max_T, max_H, max_Pi))
    else:
        with open(summary, 'a') as s:
            s.write('{},{},{},{:.2%},{:.4f},{:.4f},{:.4f},{:.6f}\n'.format(
                arg.input, rows, columns, gap_ratio, max_count, max_T, max_H,
                max_Pi))
    assert max_count > arg.resolution, (
        """
The highest resolution of given fragment is {:.2f}, which is lower than
given resolution threshold({:.2f}). Please try to use longer fragment or
lower resolution options.
""".format(max_count, arg.resolution))
    # count resolution
    (seq_count, H, Pi, T, index) = count_and_draw(alignment, consensus, arg)
    # exit if resolution lower than given threshold.
    assert len(seq_count) != 0, 'Problematic Input !'
    # find candidate
    good_region = get_good_region(index, seq_count, arg)
    consensus = find_continuous(consensus, good_region, arg.min_primer)
    primer_candidate, consensus = find_primer(consensus, arg.min_primer,
                                              arg.max_primer)
    assert len(primer_candidate) != 0, (
        'Primer not found! Try to loose options.')
    # validate
    primer_verified = validate(primer_candidate, db_file, rows, arg)
    assert len(primer_verified) != 0, 'Primer not found! Try to loose options.'
    # pick pair
    pairs = pick_pair(primer_verified, alignment, arg)
    # output
    csv_title = ('Score,Sequences,AvgProductLength,StdEV,MinProductLength,'
                 'MaxProductLength,Coverage,Resolution,TreeValue,Entropy,'
                 'LeftSeq,LeftTm,LeftAvgBitscore,LeftAvgMismatch,RightSeq,'
                 'RightTm,RightAvgBitscore,RightAvgMismatch,DeltaTm,'
                 'AlnStart,AlnEnd,AvgSeqStart,AvgSeqEnd\n')
    style = ('{:.2f},{},{:.0f},{:.0f},{},{},{:.2%},{:.2%},{:.2f},{:.2f},{},'
             '{:.2f},{:.2f},{:.2f},{},{:.2f},{:.2f},{:.2f},{:.2f},{},{},{},{}'
             '\n')
    with open('{}-{}samples-{:.2f}resolution.fastq'.format(
            arg.out, rows, arg.resolution), 'w') as out1, open(
                '{}-{}samples-{:.2f}resolution.csv'.format(
            arg.out, rows, arg.resolution), 'w') as out2:
        out2.write(csv_title)
        for pair in pairs:
            line = style.format(
                pair.score, rows, np.average(list(pair.length.values())),
                np.std(list(pair.length.values())), min(pair.length.values()),
                max(pair.length.values()), pair.coverage,
                pair.resolution, pair.tree_value, pair.entropy, pair.left.seq,
                pair.left.tm, pair.left.avg_bitscore, pair.left.avg_mismatch,
                pair.right.seq, pair.right.tm, pair.right.avg_bitscore,
                pair.right.avg_mismatch, pair.delta_tm, pair.left.start,
                pair.right.end, pair.start, pair.end)
            out2.write(line)
            SeqIO.write(pair.left, out1, 'fastq')
            SeqIO.write(pair.right, out1, 'fastq')
    print('Input alignment:')
    print('\t{}: {} rows, {} columns'.format(arg.input, rows, columns))
    print('Parameters:')
    for i in vars(arg).items():
        print('\t{}: {}'.format(i[0].capitalize(), i[1]))
    print('Found {} pairs of primers.'.format(len(pairs)))
    end = timer()
    print('Cost {:.3f} seconds.'.format(end-start))


if __name__ == '__main__':
    main()
