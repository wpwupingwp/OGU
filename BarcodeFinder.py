#!/usr/bin/python3

import argparse
import json
import re
from collections import defaultdict
from datetime import datetime
from glob import glob
from os import (devnull, environ, mkdir, pathsep, remove, sched_getaffinity,
                sep)
from os.path import abspath, basename, exists, splitext
from os.path import join as join_path
from platform import system
from shutil import unpack_archive, ReadError
from subprocess import run
from urllib import request

import numpy as np
import primer3
from Bio import Entrez, Phylo, SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline as Blast
from Bio.Data.IUPACData import ambiguous_dna_values as ambiguous_data
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from matplotlib import pyplot as plt
from matplotlib import rcParams

rcParams['lines.linewidth'] = 1.5
rcParams['axes.linewidth'] = 1.5
rcParams['axes.labelsize'] = 16
rcParams['axes.titlesize'] = 25
rcParams['font.size'] = 16


class PrimerWithInfo(SeqRecord):
    def __init__(self, seq='', quality=None, start=0, coverage=0,
                 avg_bitscore=0, mid_loc=None, avg_mismatch=0, detail=0,
                 reverse_complement=False):
        # store str
        super().__init__(Seq(seq.upper()))
        self.sequence = str(self.seq)

        """
        primer3.setGlobals seems have no effect on calcTm, so I have to replace
        all ambiguous base to G to get an approximate value. Otherwise calcTm()
        will generate -99999 if there is ambiguous base.
        """
        self.g_seq = re.sub(r'[^ATCG]', 'G', self.sequence)
        self.quality = self.letter_annotations['solexa_quality'] = quality
        self.start = self.annotations['start'] = start
        self.end = self.annotations['end'] = start + self.__len__() - 1
        self.tm = self.annotations['tm'] = primer3.calcTm(self.g_seq)
        self.coverage = self.annotations['coverage'] = coverage
        self.avg_bitscore = self.annotations['avg_bitscore'] = avg_bitscore
        self.mid_loc = self.annotations['mid_loc'] = mid_loc
        self.avg_mid_loc = 0
        self.avg_mismatch = self.annotations['avg_mismatch'] = avg_mismatch
        self.detail = self.annotations['detail'] = detail
        self.is_reverse_complement = reverse_complement
        self.description = ''
        self.hairpin_tm = 0
        self.homodimer_tm = 0
        self.update_id()

    def __getitem__(self, i):
        if isinstance(i, int):
            i = slice(i, i + 1)
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

    def reverse_complement(self, **kwargs):
        table = str.maketrans('ACGTMRWSYKVHDBXN', 'TGCAKYWSRMBDHVXN')
        new_seq = str.translate(self.sequence, table)[::-1]
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
            self.avg_mid_loc = average(list(self.mid_loc.values()))
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
        self.start = int(self.left.avg_mid_loc)
        self.end = int(self.right.avg_mid_loc)
        self.have_heterodimer = False
        self.heterodimer_tm = 0.0
        self.resolution = 0.0
        self.tree_value = 0.0
        self.entropy = 0.0
        self.pi = 0.0
        self.score = 0.0
        self.get_score()

    def __str__(self):
        return (
            'Pair(score={:.2f}, product={:.0f}, start={}, end={}, left={}, '
            'right={}, resolution={:.2%}, coverage={:.2%}, delta_tm={:.2f}, '
            'have_heterodimer={})'.format(
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
         self.tree_value) = get_resolution(alignment, self.left.start,
                                           self.right.end + 1)
        self.heterodimer_tm = primer3.calcHeterodimer(self.left.g_seq,
                                                      self.right.g_seq).tm
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
    # to be continue
    arg.add_argument('-continue', action='store_true',
                     help='continue broken download process')
    arg.add_argument('-email', default='',
                     help='email address used by NCBI Genbank')
    arg.add_argument('-fasta', default='',
                     help='unaligned fasta format data to add')
    arg.add_argument('-aln', default=None,
                     help='aligned fasta files to analyze')
    arg.add_argument('-query', help='query text')
    arg.add_argument('-stop', type=int, choices=(1, 2, 3), default=3,
                     help=('Stop after which step:\n'
                           '\t1. Download\n'
                           '\t2. Preprocess data\n'
                           '\t3. Analyze\n'))
    arg.add_argument('-no_uniq', action='store_true',
                     help='do not remove redundant records')
    output = arg.add_argument_group('Output')
    output.add_argument('-no_frag', action='store_true',
                        help='analyze whole sequence rather than divided'
                        ' fragment')
    output.add_argument('-out', help='output directory')
    output.add_argument('-rename', action='store_true',
                        help='try to rename gene')
    output.add_argument('-no_expand', default=False,
                        action='store_true',
                        help='do not expand upstream/downstream')
    output.add_argument('-expand_n', type=int, default=200,
                        help='expand length')
    filters = arg.add_argument_group('Filters')
    filters.add_argument('-gene', type=str, help='gene name')
    filters.add_argument('-group',
                         choices=('animals', 'plants', 'fungi', 'protists',
                                  'bacteria', 'archaea', 'viruses'),
                         help='Species kind')
    filters.add_argument('-min_len', default=100, type=int,
                         help='minium length')
    filters.add_argument('-max_len', default=10000, type=int,
                         help='maximum length')
    filters.add_argument('-molecular', choices=('DNA', 'RNA'),
                         help='molecular type')
    filters.add_argument('-taxon', help='Taxonomy name')
    filters.add_argument('-organelle',
                         choices=('mitochondrion', 'plastid', 'chloroplast'),
                         help='organelle type')
    options = arg.add_argument_group('Analyze')
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
    parsed = arg.parse_args()
    parsed.db_file = 'interleaved.fasta'
    parsed.no_gap_file = 'no_gap.fasta'
    if parsed.organelle is not None:
        # 10k to 1m seems enough
        parsed.min_len = 10000
        parsed.max_len = 1000000
    if (parsed.query is None and parsed.taxon is None and parsed.fasta == ''
            and parsed.aln is None and parsed.gene is None):
        arg.print_help()
        raise ValueError('Please check your input!')
    if parsed.out is None:
        parsed.out = datetime.now().isoformat().replace(':', '-')

    # arg.print_help()
    # overwrite options by given json
    if parsed.json is not None:
        with open(parsed.json, 'r') as _:
            config = json.load(_)
        d_parsed = vars(parsed)
        new_arg = argparse.Namespace(**config, **d_parsed)
        return new_arg
    else:
        return parsed


def tprint(string):
    s = '{}\t{}'.format(datetime.now().time(), string)
    print(s, flush=True)
    log_handle.write(s + '\n')


def check_tools():
    if exists('PATH.json'):
        with open('PATH.json', 'r') as path_file:
            exists_path = path_file.read().strip()
            environ['PATH'] = pathsep.join([exists_path, environ['PATH']])
    f = open(devnull, 'w')
    ok = True
    for tools in ('mafft', 'iqtree', 'blastn'):
        check = run('{} --help'.format(tools), shell=True, stdout=f, stderr=f)
        # mafft return 1 if "--help", BLAST do not have --help
        # to simplify code, use "--help" and accept 1 as returncode
        if check.returncode not in (0, 1):
            ok = False
    if not ok:
        deploy()
    f.close()


def download_software(sys):
    # althrough dirty, but save one file to make folder clean
    blast_url = ('ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/'
                 'ncbi-blast-2.7.1+')
    iqtree_url = ('https://github.com/Cibiv/IQ-TREE/release/download/v1.6.7'
                  'iqtree-1.6.7')
    mafft_url = 'https://mafft.cbrc.jp/alignment/software/mafft'
    urls = {'Linux': {'blast': blast_url+'-x64-linux.tar.gz',
                      'iqtree': iqtree_url+'-Linux.tar.gz',
                      'mafft': mafft_url+'-7.407-linux.tgz'},
            'macos': {'blast': blast_url+'.dmg',
                      'iqtree': iqtree_url+'-MacOSX.zip',
                      'mafft': mafft_url+'-7.407-mac.zip'},
            'Windows': {'blast': blast_url+'-win64.exe',
                        'iqtree': iqtree_url+'-Windows.zip',
                        'mafft': mafft_url+'-7.409-win64-signed.zip'}}
    for software in urls[sys]:
        url = urls[sys][software]
    filename = url.split('/')[-1]
    down = request.urlopen(url)
    if down.status == 200:
        with open(filename, 'wb') as out:
            out.write(down.read())
    else:
        raise Exception('Cannot download {}.'.format(software))
    try:
        unpack_archive(filename)
    except ReadError:
        pass
    return True


def deploy(software):
    sys = system()
    if sys == 'Windows':
        download_software(sys, software)
        run('ncbi-blast-2.7.1+-win64.exe', shell=True)
        environ['PATH'] = pathsep.join([
            abspath('mafft-win'), abspath('iqtree-1.6.7-Windows'+sep+'bin'),
            environ['PATH']])
    elif sys == 'Linux':
        ok = False
        for pack_mgr in ('apt', 'dnf', 'yum', 'pkg'):
            r = run('{} install ncbi-blast+ iqtree mafft'.format(pack_mgr),
                    shell=True)
            if r.returncode == 0:
                ok = True
                break
        if not ok:
            download_software(sys, software)
            environ['PATH'] = pathsep.join([
                abspath('mafft-linux64'),
                abspath('iqtree-1.6.7-Linux'+sep+'bin'),
                abspath('ncbi-blast-2.7.1+'+sep+'bin'), environ['PATH']])
    elif sys == 'Apple':
        r = run('brew --help', shell=True)
        if r.returncode != 0:
            r2 = run('/usr/bin/ruby -e "$(curl -fsSL https://raw.'
                     'githubusercontent.com/Homebrew/install/master/install)"',
                     shell=True)
            r3 = run('brew install blast mafft brewsci/science/iqtree',
                     shell=True)
            if r2.returncode != r3.returncode != 0:
                raise Exception('Cannot install brew')
        else:
            r3 = run('brew install blast mafft brewsci/science/iqtree',
                     shell=True)
        if r3.returncode != 0:
            download_software(sys)
            environ['PATH'] = pathsep.join([
                abspath('mafft-mac'), abspath('iqtree-1.6.7-MacOSX'+sep+'bin'),
                abspath('ncbi-blast-2.7.1+'+sep+'bin'), environ['PATH']])
    with open('PATH.json', 'w') as path_out:
        json.dump(environ['PATH'], path_out)
    return environ['PATH']



def get_query_string(arg):
    condition = list()
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
    tprint('Your query:\n\t{}'.format(query))
    Entrez.email = arg.email
    query_handle = Entrez.read(Entrez.esearch(db='nuccore', term=query,
                                              usehistory='y'))
    count = int(query_handle['Count'])
    tprint('{} records found.'.format(count))
    tprint('Downloading... Ctrl+C to quit')
    json_file = join_path(arg.out, 'Query.json')
    with open(json_file, 'w') as _:
        json.dump(query_handle, _, indent=4, sort_keys=True)
    if arg.query is None:
        name = safe(arg.taxon)
    else:
        name = safe(arg.query)
    file_name = join_path(arg.out, name + '.gb')
    output = open(file_name, 'w')
    ret_start = 0
    ret_max = 1000
    while ret_start <= count:
        tprint('{:4d}-{:4d}'.format(ret_start, ret_start + ret_max))
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
        ret_start += 1000
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
        try:
            codon = Seq(re.search(pattern, lower).group(1))
        except AttributeError:
            return old_name, 'bad_name'
        try:
            new_name = 'trn{}{}'.format(codon.reverse_complement().translate(),
                                        codon.transcribe())
        except ValueError:
            return old_name, 'bad_name'
        gene_type = 'tRNA'
    elif lower.startswith('rrn'):
        pattern = re.compile(r'(\d+\.?\d?)')
        try:
            number = re.search(pattern, lower).group(1)
        except AttributeError:
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
        try:
            gene = match.group('gene')
            suffix = match.group('suffix')
        except ValueError:
            return old_name, 'bad_name'
        new_name = '{}{}'.format(gene, suffix.upper())
        # captitalize last letter
        if len(new_name) > 3:
            s = list(new_name)
            if s[-1].isalpha():
                new_name = '{}{}'.format(
                    ''.join(s[:-1]), ''.join(s[-1]).upper())
        gene_type = 'normal'
    if len(lower) >= 15:
        gene_type = 'suspicious_name'
    return new_name, gene_type


def safe(old):
    return re.sub(r'\W', '_', old)


def get_taxon(order_family):
    """
    From Zhang guojin
    order end with ales
    family end with aceae except 8
    http://duocet.ibiodiversity.net/index.php?title=%E4%BA%92%E7%94%A8%E5%90%8D
    %E7%A7%B0&mobileaction=toggle_view_mobile"""
    # order|family|organims(genus|species)
    family_exception_raw = (
        'Umbelliferae,Palmae,Compositae,Cruciferae,Guttiferae,Leguminosae,'
        'Leguminosae,Papilionaceae,Labiatae,Gramineae')
    family_exception = family_exception_raw[0].split(',')
    order = ''
    family = ''
    for item in order_family:
        if item.endswith('ales'):
            order = item
        elif item.endswith('aceae'):
            family = item
        elif item in family_exception:
            family = item
    return order, family


def write_seq(name, sequence_id, feature, whole_seq, path, arg):
    """
    Write fasta file.
    """
    filename = join_path(path, name + '.fasta')
    sequence = feature.extract(whole_seq)

    with open(filename, 'a') as handle:
        handle.write(sequence_id + '\n')
        handle.write(str(sequence) + '\n')
    if not arg.no_expand:
        if feature.location_operator == 'join':
            loc = feature.location.parts
            # ensure increasing order
            # parts do not have sort method
            sorted(loc, key=lambda x: x.start)
            new_loc = sum([
                # avoid IndexError
                FeatureLocation(max(0, loc[0].start - arg.expand_n),
                                loc[0].end, loc[0].strand),
                *loc[1:-1],
                FeatureLocation(loc[-1].start,
                                min(len(whole_seq), loc[-1].end+arg.expand_n),
                                loc[-1].strand)])
            feature.location = new_loc
        feature.type = 'expand'
        sequence = feature.extract(whole_seq)
        filename2 = join_path(path, 'expand.{}.fasta'.format(name))
        with open(filename2, 'a') as handle:
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
                                'intergenic_spacer')
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
    spacers = list()
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
    groupby_gene = '{}-groupby_gene'.format(gbfile.replace('.gb', ''))
    mkdir(groupby_gene)
    groupby_name = '{}-groupby_name'.format(gbfile.replace('.gb', ''))
    mkdir(groupby_name)
    handle_raw = open(gbfile + '.fasta', 'w')
    wrote_by_gene = set()
    wrote_by_name = set()

    for record in SeqIO.parse(gbfile, 'gb'):
        # only accept gene, product, and spacer in misc_features.note
        order_family = record.annotations['taxonomy']
        order, family = get_taxon(order_family)
        organism = record.annotations['organism'].replace(' ', '_')
        genus, *species = organism.split('_')
        taxon = '{}|{}|{}|{}'.format(order, family, genus, '_'.join(species))
        accession = record.annotations['accessions'][0]
        try:
            specimen = record.features[0].qualifiers['specimen_voucher'
                                                     ][0].replace(' ', '_')
        except (IndexError, KeyError):
            specimen = ''
        whole_seq = record.seq
        feature_name = list()
        genes = list()

        for feature in record.features:
            name, feature_type = get_feature_name(feature, arg)
            # skip unsupport feature
            if name is None:
                continue
            if len(name) > 50:
                tprint('Too long name: {}.'.format(name))
                name = name[:50] + '...'
            # skip abnormal annotation
            if len(feature) > 20000:
                tprint('Skip abnormal annotaion of {}(Accession {})!'.format(
                    name, accession))
                continue
            if feature_type == 'gene':
                genes.append([name, feature])
            feature_name.append(name)
            sequence_id = '>' + '|'.join([name, taxon, accession, specimen])
            wrote = write_seq(name, sequence_id, feature, whole_seq,
                              groupby_gene, arg)
            wrote_by_gene.add(wrote)

        # extract spacer
        spacers = get_spacer(genes)
        for spacer in spacers:
            sequence_id = '>' + '|'.join([spacer.id, taxon,
                                          accession, specimen])
            wrote = write_seq(spacer.id, sequence_id, spacer, whole_seq,
                              groupby_gene, arg)
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
        filename = join_path(groupby_name, name_str + '.fasta')
        with open(filename, 'a') as out:
            SeqIO.write(record, out, 'fasta')
            wrote_by_name.add(filename)
        # write raw fasta
        SeqIO.write(record, handle_raw, 'fasta')

    # skip analyze of Unknown.fasta
    unknown = join_path(groupby_name, 'Unknown.fasta')
    if unknown in wrote_by_name:
        wrote_by_name.remove(unknown)
    tprint('Divide done.')
    return list(wrote_by_gene), list(wrote_by_name)


def uniq(files):
    """
    Only keep first record with same id.
    """
    uniq_files = list()
    for fasta in files:
        raw = SeqIO.parse(fasta, 'fasta')
        new = fasta + '.uniq'
        new_handle = open(new, 'w')
        names = set()
        for record in raw:
            # gene|order|family|genus|species|specimen
            name = ' '.join(record.id.split('|')[3:5])
            if name in names:
                pass
            else:
                SeqIO.write(record, new_handle, 'fasta')
                names.add(name)
        new_handle.close()
        uniq_files.append(new)
    return uniq_files


def mafft(files):
    result = list()
    # get available CPU cores
    cores = len(sched_getaffinity(0))
    for fasta in files:
        tprint('Aligning {}'.format(fasta))
        out = fasta + '.aln'
        _ = ('mafft --thread {} --reorder --quiet --adjustdirection '
             '{} > {}'.format(cores - 1, fasta, out))
        m = run(_, shell=True)
        if m.returncode == 0:
            result.append(out)
    tprint('Alignment done.')
    for i in glob('_order*'):
        remove(i)
    return result


def average(x):
    if len(x) == 0:
        return 0
    else:
        return sum(x) / len(x)


def prepare(arg):
    """
    Given fasta format alignment filename, return a numpy array for sequence:
    Generate fasta file without gap for makeblastdb, return file name.
    """
    data = list()
    record = ['id', 'sequence']
    with open(arg.input, 'r') as raw, open(arg.no_gap_file, 'w') as no_gap:
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
        tprint('{} does not have uniform width!'.format(arg.input))
        return None, None, None

    # Convert List to numpy array.
    # order 'F' is a bit faster than 'C'
    # new = np.hstack((name, seq)) -> is slower
    name = np.array([[i[0]] for i in data], dtype=np.bytes_)
    sequence = np.array([list(i[1]) for i in data], dtype=np.bytes_, order='F')

    if name is None:
        tprint('Bad fasta file {}.'.format(arg.input))
        name = None
    # tree require more than 4 sequences
    if len(sequence) < 4:
        tprint('Too few sequence in {} (less than 4)!'.format(arg.input))
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
    frequency = list()
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
        with open(devnull, 'w') as f:
            iqtree = run('iqtree -s {} -m JC -fast -czb'.format(aln_file),
                         stdout=f, stderr=f, shell=True)
        # just return 0 if there is error
        if iqtree.returncode != 0:
            tprint('Cannot get tree_value of region {}-{} bp!'.format(
                start, end))
            clean()
            return resolution, entropy, pi, 0
        tree = Phylo.read(aln_file + '.treefile', 'newick')
        # skip the first empty node
        internals = tree.get_nonterminals()[1:]
        tree_value = len(internals) / len(tree.get_terminals())
        clean()
    else:
        tree_value = 0
    return resolution, entropy, pi, tree_value


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
    primers = list()
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
    http://www.tiem.utk.edu/~gross/bioed/bealsmodules/shannonDI.html
    return List[float]
    All calculation excludes primer sequence.
    """
    rows, columns = alignment.shape
    min_primer = arg.min_primer
    max_product = arg.max_product
    step = arg.step
    out_file = join_path(arg.out, basename(arg.out_file))
    # r_list, h_list, pi_list, t_list : count, normalized entropy, Pi and
    #  tree value
    r_list = list()
    h_list = list()
    pi_list = list()
    t_list = list()
    max_h = np.log2(rows)
    index = list()
    max_plus = max_product - min_primer * 2
    max_range = columns - max_product
    for i in range(0, max_range, step):
        # skip gap
        # if consensus.sequence[i] in ('-', 'N'):
        #     continue
        # exclude primer sequence
        resolution, entropy, pi, tree_value = get_resolution(
            alignment, i, i + max_plus, arg.fast)
        r_list.append(resolution)
        h_list.append(entropy / max_h)
        pi_list.append(pi)
        t_list.append(tree_value)
        index.append(i)

    plt.style.use('seaborn-colorblind')
    fig, ax1 = plt.subplots(figsize=(15 + len(index) // 5000, 10))
    plt.title('Resolution(window={} bp, step={} bp)\n'.format(
        max_product, step))
    plt.xlabel('Base')
    # plt.xticks(np.linspace(0, max_range, 21))
    if not arg.fast:
        ax1.set_ylabel('Normalized Shannon Index / Resolution / TreeValue')
        ax1.plot(index, t_list, label='TreeValue', alpha=0.8)
    else:
        ax1.set_ylabel('Normalized Shannon Index / Resolution')

    ax1.plot(index, h_list, label='Shannon Index', alpha=0.8)
    ax1.plot(index, r_list, label='Resolution', alpha=0.8)
    ax1.legend(loc='lower left')
    ax1.yaxis.set_ticks(np.linspace(0, 1, num=11))
    ax2 = ax1.twinx()
    ax2.plot(index, pi_list, 'k-', label=r'$\pi$', alpha=0.8)
    ax2.set_ylabel(r'$\pi$', rotation=-90, labelpad=20)
    # _ = round(np.log10(max(Pi))) + 1
    # ax2.yaxis.set_ticks(np.linspace(0, 10**_, num=11))
    ax2.legend(loc='upper right')
    # plt.yscale('log')
    plt.savefig(out_file + '.pdf')
    plt.savefig(out_file + '.png')
    # plt.show()
    with open(out_file + '-Resolution.tsv', 'w') as _:
        _.write('Index,R,H,Pi,T\n')
        for i, r, h, pi, t in zip(index, r_list, h_list, pi_list, t_list):
            _.write('{},{:.2f},{:.2f},{:.2f},{:.2f}\n'.format(i, r, h, pi, t))
    return r_list, h_list, pi_list, t_list, index


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
    query_file = arg.out_file + '.candidate.fasta'
    # SeqIO.write fasta file directly is prohibited. have to write fastq at
    with open(query_file + '.fastq', 'w') as _:
        SeqIO.write(primer_candidate, _, 'fastq')
    SeqIO.convert(query_file + '.fastq', 'fastq', query_file, 'fasta')
    # build blast db
    with open(devnull, 'w') as f:
        _ = run('makeblastdb -in {} -dbtype nucl'.format(db_file),
                shell=True, stdout=f)
        if _.returncode != 0:
            tprint('makeblastdb error!')
            return []
    # blast
    tprint('Validate with BLAST.')
    blast_result_file = 'blast.result.tsv'
    fmt = 'qseqid sseqid qseq nident mismatch score qstart qend sstart send'
    cmd = Blast(num_threads=len(sched_getaffinity(0)),
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
    # clean makeblastdb files
    for i in glob(db_file + '*'):
        remove(i)
    remove(blast_result_file)
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
    less_pairs = list()
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
    tprint('{} pairs of redundant primers were removed'.format(
        len(pairs) - len(less_pairs)))
    good_pairs = list()
    for i in less_pairs:
        i.add_info(alignment)
        if i.resolution >= arg.resolution:
            good_pairs.append(i)
    good_pairs.sort(key=lambda x: x.score, reverse=True)
    tprint('Successfully found {} pairs of validated primers'.format(
        len(good_pairs)))
    return good_pairs


def analyze(arg):
    """
    Automatic design primer for DNA barcode.
    """
    arg.out_file = splitext(arg.input)[0]
    # read from fasta, generate new fasta for makeblastdb
    name, alignment, db_file = prepare(arg)
    if name is None:
        return
    rows, columns = alignment.shape
    # generate consensus
    base_cumulative_frequency = count_base(alignment, rows, columns)
    tprint('Generate consensus.')
    consensus = generate_consensus(base_cumulative_frequency, arg.coverage,
                                   rows, arg.out_file + '.consensus.fastq')
    tprint('Evaluate whole alignment.')
    max_count, max_h, max_pi, t = get_resolution(alignment, 0, columns)
    n_gap = sum([i[5] for i in base_cumulative_frequency])
    gap_ratio = n_gap / rows / columns
    summary = join_path(arg.out, 'Summary.csv')
    if not exists(summary):
        with open(summary, 'w') as s:
            s.write('Name,Sequences,Length,GapRatio,ObservedResolution,'
                    'TreeValue,ShannonIndex,Pi\n')
            s.write('{},{},{},{:.2%},{:.6f},{:.6f},{:.6f},{:.6f}\n'.format(
                basename(arg.input), rows, columns, gap_ratio,
                max_count, t, max_h, max_pi))
    else:
        with open(summary, 'a') as s:
            s.write('{},{},{},{:.2%},{:.4f},{:.4f},{:.4f},{:.6f}\n'.format(
                basename(arg.input), rows, columns, gap_ratio,
                max_count, t, max_h, max_pi))
    if max_count < arg.resolution:
        tprint('Too low resolution of {} !'.format(arg.input))
        return
    # count resolution
    tprint('Sliding window analyze.')
    (seq_count, H, Pi, T, index) = count_and_draw(alignment, arg)
    # exit if resolution lower than given threshold.
    if len(seq_count) == 0:
        tprint('Problematic Input of {}.!'.format(arg.input))
    # find candidate
    good_region = get_good_region(index, seq_count, arg)
    consensus = find_continuous(consensus, good_region, arg.min_primer)
    tprint('Find candidate primer pairs')
    primer_candidate, consensus = find_primer(consensus, arg.min_primer,
                                              arg.max_primer)
    if len(primer_candidate) == 0:
        tprint('Cannot find primers in {}. Try to loose options!'.format(
            arg.input))
        return
    tprint('Found {} candidate primers'.format(len(primer_candidate)))
    # validate
    primer_verified = validate(primer_candidate, db_file, rows, arg)
    if len(primer_verified) == 0:
        tprint('Cannot find primers in {}. Try to loose options!'.format(
            arg.input))
        return
    # pick pair
    pairs = pick_pair(primer_verified, alignment, arg)
    if len(pairs) == 0:
        tprint('Cannot find primers in {}. Try to loose options!'.format(
            arg.input))
        return
    # output
    csv_title = ('Score,Sequences,AvgProductLength,StdEV,MinProductLength,'
                 'MaxProductLength,Coverage,Resolution,TreeValue,Entropy,'
                 'LeftSeq,LeftTm,LeftAvgBitscore,LeftAvgMismatch,RightSeq,'
                 'RightTm,RightAvgBitscore,RightAvgMismatch,DeltaTm,'
                 'AlnStart,AlnEnd,AvgSeqStart,AvgSeqEnd\n')
    style = ('{:.2f},{},{:.0f},{:.0f},{},{},{:.2%},{:.2%},{:.2f},{:.2f},{},'
             '{:.2f},{:.2f},{:.2f},{},{:.2f},{:.2f},{:.2f},{:.2f},{},{},{},{}'
             '\n')
    _ = join_path(arg.out, basename(arg.out_file))
    with open(_ + '.fastq', 'w') as out1, open(_ + '.csv', 'w') as out2:
        out2.write(csv_title)
        for pair in pairs:
            line = style.format(
                pair.score, rows, average(list(pair.length.values())),
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
    tprint('Primers info were written into {}.csv'.format(_))


def main():
    """
    1. Get data from Genbank.
    2. Divide according to annotation.
    3. Analyze.
    """
    arg = parse_args()
    mkdir(arg.out)
    global log_handle
    log_handle = open(join_path(arg.out, 'Log.txt'), 'w')
    tprint('Welcome to BarcodeFinder!')
    check_tools()

    def analyze_wrapper(files):
        for aln in files:
            tprint('Analyze {}.'.format(aln))
            arg.input = aln
            analyze(arg)
        # dirty work
        try:
            remove(arg.no_gap_file)
            remove(arg.db_file)
        except FileNotFoundError:
            pass

    if arg.aln is not None:
        analyze_wrapper(list(glob(arg.aln)))
        return
    user_data = list(glob(arg.fasta))
    query = get_query_string(arg)
    tprint('Generate query for Genbank.')
    if query is not None:
        tprint('Download data from Genbank.')
        gbfile = download(arg, query)
        tprint('Divide data by annotation.')
        wrote_by_gene, wrote_by_name = divide(gbfile, arg)
        wrote_by_gene.extend(user_data)
        wrote_by_name.extend(user_data)
    else:
        wrote_by_gene = user_data[::]
        wrote_by_name = user_data[::]
    if arg.stop == 1:
        return
    if arg.no_uniq:
        pass
    else:
        tprint('Remove redundant sequences.')
        wrote_by_gene = uniq(wrote_by_gene)
        wrote_by_name = uniq(wrote_by_name)
    tprint('Aligning sequences.')
    if not arg.no_frag or arg.max_len > 10000:
        # less than two records will cause empty output, which was omit
        aligned = mafft(wrote_by_gene)
    else:
        aligned = mafft(wrote_by_name)
    if arg.stop == 2:
        return
    analyze_wrapper(aligned)
    tprint('Finished. You can find output in {}.'.format(arg.out))
    tprint('Summary info were written into {}.'.format(
        join_path(arg.out, 'Summary.csv')))
    _ = join_path(arg.out, 'Options.json')
    with open(_, 'w') as out:
        json.dump(vars(arg), out, indent=4, sort_keys=True)
    tprint('Options were dumped into {}.'.format(_))
    log_handle.close()
    return


if __name__ == '__main__':
    main()
