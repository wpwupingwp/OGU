#!/usr/bin/python3

import argparse
import logging
import re
from collections import defaultdict
from itertools import product as cartesian_product
from os import cpu_count, devnull
from pathlib import Path
from subprocess import run

import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Data.IUPACData import ambiguous_dna_values as ambiguous_data
from Bio.Blast.Applications import NcbiblastnCommandline as Blast
from Bio import SeqIO
from primer3 import (calcTm, calcHairpinTm, calcHomodimerTm,
                     calcHeterodimerTm)

from BarcodeFinder import utils
from BarcodeFinder import evaluate

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


class Pair:
    # save memory
    __slots__ = ['left', 'right', 'delta_tm', 'coverage', 'start', 'end',
                 'observed_res', 'tree_res', 'pd_terminal', 'entropy',
                 'have_heterodimer', 'heterodimer_tm', 'pi', 'score',
                 'length', 'gap_ratio']


    _title = ('Score,AvgProductLength,StdEV,'
              'MinProductLength,MaxProductLength,'
              'Coverage,Observed_Res,Tree_Res,PD_terminal,Entropy,'
              'LeftSeq,LeftTm,LeftAvgBitscore,LeftAvgMismatch,'
              'RightSeq,RightTm,RightAvgBitscore,RightAvgMismatch,'
              'DeltaTm,AlnStart,AlnEnd,AvgSeqStart,AvgSeqEnd')
    _style = ('{:.2f},{:.0f},{:.0f},'
              '{},{},'
              '{:.2%},{:.2%},{:.6f},{:.6f},{:.6f},'
              '{},{:.2f},{:.2f},{:.2f},'
              '{},{:.2f},{:.2f},{:.2f},'
              '{:.2f},{},{},{},{}')

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
        self.observed_res = 0.0
        self.tree_res = 0.0
        self.pd_terminal = 0.0
        self.entropy = 0.0
        self.pi = 0.0
        self.gap_ratio = 0.0
        self.score = self.get_score()

    def __repr__(self):
        return (
            'Pair(score={:.2f}, product={:.0f}, start={}, end={}, left={}, '
            'right={}, observerd_res={:.2%}, coverage={:.2%},'
            'delta_tm={:.2f}, have_heterodimer={})'.format(
                self.score, utils.safe_average(
                    list(self.length.values())), self.start,
                self.end, self.left.seq, self.right.seq, self.observed_res,
                self.coverage, self.delta_tm, self.have_heterodimer))

    def __str__(self):
        return Pair._style.format(
            self.score, utils.safe_average(list(self.length.values())),
            np.std(list(self.length.values())), min(self.length.values()),
            max(self.length.values()),
            self.coverage, self.observed_res, self.tree_res,
            self.pd_terminal, self.entropy,
            self.left.seq, self.left.tm, self.left.avg_bitscore,
            self.left.avg_mismatch,
            self.right.seq, self.right.tm, self.right.avg_bitscore,
            self.right.avg_mismatch,
            self.delta_tm, self.left.start, self.right.end, self.start,
            self.end)

    def get_score(self):
            # calculate score of given primer pairs. Suggestion only
        # use score to filter primer pairs
        return (utils.safe_average(list(self.length.values())) * 0.5
                + self.coverage * 200
                + len(self.left) * 10
                + len(self.right) * 10
                + self.observed_res * 100
                + self.tree_res * 100 + self.entropy * 5
                - int(self.have_heterodimer) * 10
                - self.delta_tm * 5 - self.left.avg_mismatch * 10
                - self.right.avg_mismatch * 10)

    def add_info(self, alignment):
        # put attributes that need heavy computation here for the final primer
        # pairs in order to save CPU time
        if not self.right.is_reverse_complement:
            self.right = self.right.reverse_complement()
        # include end base, use alignment loc for slice
        subalign = alignment[:, self.left.start:self.right.end+1]
        tmp = Path()
        variance, _ = evaluate.get_resolution(subalign, tmp)
        (self.gap_ratio, self.observed_res, self.entropy, self.pi,
         _, _, self.pd_terminal) = variance[2:9]
        self.heterodimer_tm = calc_ambiguous_seq(
            calcHeterodimerTm, self.left.seq, self.right.seq)
        if max(self.heterodimer_tm, self.left.tm,
               self.right.tm) == self.heterodimer_tm:
            self.have_heterodimer = True
        else:
            self.have_heterodimer = False
        self.get_score()
        self.left.update_id()
        self.right.update_id()
        return self


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


def parse_args(arg_str=None):
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=primer_main.__doc__)
    arg.add_argument('-aln', nargs='*', help='alignment files')
    arg.add_argument('-aln_folder', default=None,
                     help='folder of aligned files')
    arg.add_argument('-out', help='output directory')
    arg.add_argument('-a', dest='ambiguous_base_n', default=4, type=int,
                        help='number of ambiguous bases')
    arg.add_argument('-c', dest='coverage', default=0.6, type=float,
                        help='minimal coverage of base and primer')
    arg.add_argument('-m', dest='mismatch', default=4, type=int,
                        help='maximum mismatch bases in primer')
    arg.add_argument('-pmin', dest='min_primer', default=20, type=int,
                        help='minimum primer length')
    arg.add_argument('-pmax', dest='max_primer', default=25, type=int,
                        help='maximum primer length')
    arg.add_argument('-r', dest='resolution', type=float, default=0.5,
                        help='minimal resolution')
    arg.add_argument('-t', dest='top_n', type=int, default=1,
                        help='keep n primers for each high variant region')
    arg.add_argument('-tmin', dest='min_product', default=350, type=int,
                        help='minimum product length(include primer)')
    arg.add_argument('-tmax', dest='max_product', default=600, type=int,
                        help='maximum product length(include primer)')
    if arg_str is None:
        return arg.parse_args()
    else:
        return arg.parse_known_args(arg_str.split(' '))[0]


def init_arg(arg):
    if arg.aln is None and arg.aln_folder is None:
        log.error('Empty input.')
        return None
    if all([arg.aln, arg.aln_folder]):
        log.critical('Cannot use "-aln" and "-aln_folder" at same time!')
        log.critical('Ignore "-aln" option.')
        arg.aln = None
    if arg.aln is not None:
        arg.aln = [Path(i).absolute() for i in arg.aln]
    if arg.fasta_folder is not None:
        # overwrite
        arg.aln = [i.absolute() for i in Path(arg.aln_folder).glob('*')]
    arg = utils.init_out(arg)
    if arg.out is None:
        return None
    return arg


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
        log.warning('Too many ambiguous bases. Skip')
        return 0
    seq_str = _expand(seq)
    if seq2 is None:
        values = [func(i) for i in seq_str]
    else:
        if len(seq2) > LEN_LIMIT:
            log.warning('Too many ambiguous bases. Skip')
            return 0
        seq_str2 = _expand(seq2)
        products = cartesian_product(seq_str, seq_str2)
        values = [func(i[0], i[1]) for i in products]
    # primer3 will return negative values sometime
    values_positive = [max(0, i) for i in values]
    return utils.safe_average(values_positive)


def count_base(alignment: np.array):
    """
    Given alignment numpy array, count cumulative frequency of base in each
    column (consider ambiguous base and "N", "-" and "?", otherwise omit).
    Return [[float, float, float, float, float, float, float]] for
    [A, T, C, G, N, GAP, OTHER].
    """
    frequency = []
    rows, columns = alignment.shape
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


def get_consensus(base_cumulative_frequency, coverage_percent: float,
                  rows: int, output: Path):
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


def validate(primer_candidate, locus_name, n_seqs, arg):
    """
    Do BLAST. Parse BLAST result. Return list of PrimerWithInfo which passed
    the validation.
    """
    EVALUE = 1e-2
    query_file = arg._primer / (locus_name+'.candidate.fasta')
    query_file_fastq = arg._primer / (locus_name+'.candidate.fastq')
    # SeqIO.write fasta file directly is prohibited. have to write fastq at
    # first.
    with open(query_file_fastq, 'w', encoding='utf-8') as _:
        SeqIO.write(primer_candidate, _, 'fastq')
    SeqIO.convert(query_file_fastq, 'fastq', query_file, 'fasta')
    # build blast db
    db_file = arg._tmp / (locus_name+'db_file.fasta')
    with open(devnull, 'w', encoding='utf-8') as f:
        _ = run(f'makeblastdb -in {db_file} -dbtype nucl', shell=True, stdout=f)
        if _.returncode != 0:
            log.critical('Failed to run makeblastdb. Skip BLAST.')
            return []
    # BLAST
    blast_result_file = arg._tmp / (locus_name+'blast.result.tsv')
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
    for query in utils.parse_blast_tab(blast_result_file):
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
            loc = utils.safe_average([hit.hit_start, hit.hit_end])
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
    # clean
    utils.clean_tmp(db_file)
    blast_result_file.unlink()
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
            cluster = []
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


def get_observed_res(alignment: np.array, arg):
    """
    For primer design.
    Observed resolution used as lower bound.
    Args:
        alignment:
        arg:
    Returns:
        index: list
        observed_res_list: list
    """
    index = []
    observed_res_list = []
    rows, columns = alignment.shape
    for i in range(0, columns, arg.step):
        subalign = alignment[:, i+arg.size]
        sub_rows, sub_columns = subalign.shape
        item, count = np.unique(subalign, return_counts=True, axis=0)
        observed_res = len(count) / sub_rows
        index.append(i)
        observed_res_list.append(observed_res)
    return index, observed_res_list


def primer_design(aln: Path, result: Path, arg):
    locus_name = aln.stem
    name, alignment, = evaluate.fasta_to_array(aln)
    if name is None:
        log.info(f'Invalid alignment file {aln}.')
        return False
    rows, columns = alignment.shape
    # generate consensus
    base_cumulative_frequency = count_base(alignment)
    log.info(f'Generate consensus of {aln}.')
    consensus = get_consensus(base_cumulative_frequency, arg.coverage, rows,
                              arg.out/(locus_name+'.consensus.fastq'))
    log.info('Evaluate whole alignment of {aln}.')
    # a_ : alignment
    item, count = np.unique(alignment, return_counts=True, axis=0)
    observed_res = len(count) / rows
    if observed_res < arg.resolution:
        log.error('Observed resolution is too low.')
        return False
    log.info('Start the sliding-window scan.')
    observed_res_list, index = get_observed_res(alignment, arg)
    log.info('Evaluation finished.')
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
    log.info(f'Found {len(primer_candidate)} candidate primers.')
    log.info('Validate with BLAST. May be slow.')
    primer_verified = validate(primer_candidate, locus_name, rows, arg)
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
    out1 = open(arg._primer/(locus_name+'.primer.fastq'), 'w', encoding='utf-8')
    out2 = open(arg._primer/(locus_name+'.primer.csv'), 'w', encoding='utf-8')
    # write primers to one file
    out3 = open(result, 'a', encoding='utf-8')
    csv_title = 'Locus,Samples,' + Pair._title + '\n'
    out2.write(csv_title)
    for pair in pairs:
        line = f'{locus_name},{rows},{str(pair)}\n'
        out2.write(line)
        out3.write(line)
        SeqIO.write(pair.left, out1, 'fastq')
        SeqIO.write(pair.right, out1, 'fastq')
    out1.close()
    out2.close()
    out3.close()
    log.info(f'Primers info were written into {out2}')
    return True


def primer_main(arg_str=None):
    """
    Evaluate variance of alignments.
    Args:
        arg_str:
    Returns:
        aln: aligned files
        out_csv: evaluation of each locus
    """
    log.info('Running primer module...')
    arg = parse_args(arg_str)
    arg = init_arg(arg)
    if arg is None:
        log.error('Quit.')
        return None
    primer_result = arg.out / 'Primers.csv'
    csv_title = 'Locus,Samples,' + Pair._title + '\n'
    with open(primer_result, 'w', encoding='utf-8') as out:
        out.write(csv_title)
    if arg is None:
        log.info('Quit.')
        return None
    for aln in arg.aln:
        primer_design(aln, primer_result, arg)

    log.info(f'Primer result could be found in {primer_result}')
    log.info('Primer module Finished.')
    return arg, arg.aln


if __name__ == '__main__':
    primer_main()