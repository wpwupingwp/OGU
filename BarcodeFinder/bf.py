#!/usr/bin/python3

import argparse
import logging

from BarcodeFinder import utils
from BarcodeFinder import gb2fasta
from BarcodeFinder import evaluate
from BarcodeFinder import primer


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
        description=bf_main.__doc__)
    general = arg.add_argument_group('General')
    general.add_argument('-aln', help='aligned fasta files to analyze')
    general.add_argument('-fasta', help='unaligned fasta format data to add')
    general.add_argument('-gb', help='genbank files')
    general.add_argument('-out', help='output directory')
    gb2fasta_ = arg.add_argument_group('GB2Fasta')
    # genes in IR regions
    gb2fasta_.add_argument('-allow_mosaic_spacer', action='store_true',
                           help='allow mosaic spacer')
    gb2fasta_.add_argument('-allow_repeat', action='store_true',
                           help='allow repeat genes or spacer')
    gb2fasta_.add_argument('-allow_invert_repeat', action='store_true',
                           help='allow invert-repeat spacers')
    # for primer
    gb2fasta_.add_argument('-expand', type=int, default=0,
                           help='expand length of upstream/downstream')
    gb2fasta_.add_argument('-max_name_len', default=100, type=int,
                           help='maximum length of feature name')
    # handle rps12
    gb2fasta_.add_argument('-max_seq_len', default=20000, type=int,
                           help='maximum length of feature sequence')
    gb2fasta_.add_argument('-no_divide', action='store_true',
                           help='only download')
    # for plastid genes
    gb2fasta_.add_argument('-rename', action='store_true',
                           help='try to rename gene')
    gb2fasta_.add_argument('-unique', choices=('longest', 'first', 'no'),
                           default='first',
                           help='method to remove redundant sequences')
    gb2fasta_.add_argument('-email', type=str,
                           help='email address for querying Genbank')
    gb2fasta_.add_argument('-exclude', type=str, help='exclude option')
    gb2fasta_.add_argument('-gene', type=str, help='gene name')
    # in case of same taxonomy name in different group
    gb2fasta_.add_argument('-group',
                           choices=('all', 'animals', 'plants', 'fungi',
                                    'protists', 'bacteria', 'archaea',
                                    'viruses'),
                           default='all',
                           help='Species kind')
    gb2fasta_.add_argument('-min_len', default=100, type=int,
                           help='minimum length')
    gb2fasta_.add_argument('-max_len', default=10000, type=int,
                           help='maximum length')
    gb2fasta_.add_argument('-date_start', type=str,
                           help='release date beginning, (eg. 1970/1/1)')
    gb2fasta_.add_argument('-date_end', type=str,
                           help='release date end, (eg. 2020/12/31)')
    gb2fasta_.add_argument('-molecular', choices=('all', 'DNA', 'RNA'),
                           default='all', help='molecular type')
    gb2fasta_.add_argument('-og', '-organelle', dest='organelle',
                           choices=('ignore', 'both', 'no', 'mt',
                                    'mitochondrion', 'cp', 'chloroplast',
                                    'pl', 'plastid'),
                           default='ignore', help='organelle type')
    gb2fasta_.add_argument('-query', nargs='*', help='query text')
    gb2fasta_.add_argument('-refseq', choices=('both', 'yes', 'no'),
                           default='both', help='include RefSeq or not')
    gb2fasta_.add_argument('-seq_n', default=0, type=int,
                           help='maximum number of records to download')
    gb2fasta_.add_argument('-taxon', help='Taxonomy name')
    evaluate = arg.add_argument_group('Evaluate')
    evaluate.add_argument('-ig', '-ignore_gap', dest='ignore_gap',
                          action='store_true',
                          help='ignore gaps in alignments')
    evaluate.add_argument('-iab', '-ignore_ambiguous',
                          dest='ignore_ambiguous_base', action='store_true',
                          help='ignore ambiguous bases like "M" or "N"')
    evaluate.add_argument('-quick', action='store_true',
                          help='skip sliding-window analysis')
    evaluate.add_argument('-size', type=int, default=500,
                          help='window size')
    evaluate.add_argument('-step', default=50, type=int,
                          help='step length for sliding-window scan')
    evaluate.add_argument('-skip_primer', action='store_true',
                          help='skip primer designing')
    primer = arg.add_argument_group('Primer')
    primer.add_argument('-ambiguous', dest='ambiguous_base_n', default=4,
                        type=int, help='number of ambiguous bases')
    primer.add_argument('-coverage', dest='coverage', default=0.5, type=float,
                        help='minimal coverage of base and primer')
    primer.add_argument('-mismatch', dest='mismatch', default=4, type=int,
                        help='maximum mismatch bases in primer')
    primer.add_argument('-pmin', dest='min_primer', default=20, type=int,
                        help='minimum primer length')
    primer.add_argument('-pmax', dest='max_primer', default=25, type=int,
                        help='maximum primer length')
    primer.add_argument('-res', dest='resolution', type=float, default=0.3,
                        help='minimal resolution')
    primer.add_argument('-topn', dest='top_n', type=int, default=1,
                        help='keep n primers for each high variant region')
    primer.add_argument('-tmin', dest='min_product', default=350, type=int,
                        help='minimum product length(include primer)')
    primer.add_argument('-tmax', dest='max_product', default=600, type=int,
                        help='maximum product length(include primer)')
    return arg.parse_args()


def init_arg(arg):
    utils.get_all_third_party()
    arg = utils.init_out(arg, from_main=True)
    if arg.out is None:
        return None
    query = gb2fasta.get_query_string(arg, silence=True)
    if not any([arg.gb, arg.fasta, arg.aln, arg.query, query]):
        log.error('Empty input.')
        return None
    return arg


def bf_main():
    """
    Call gb2fasta, evaluate and primer functions.
    """
    log.info('Welcome to BarcodeFinder.')
    arg = parse_args()
    arg = init_arg(arg)
    if arg is None:
        log.error('Quit.')
        return
    log_file = arg.out / 'Log.txt'
    log_file_handler = logging.FileHandler(log_file, mode='a')
    log_file_handler.setLevel(logging.INFO)
    log.addHandler(log_file_handler)

    option = utils.arg_to_str(arg)
    arg, other_args, = gb2fasta.gb2fasta_main()
    arg, other_args2 = evaluate.evaluate_main()
    arg, other_args3 = primer.primer_main()
    if len(other_args3) != 0:
        log.debug(f'Unrecognized options {" ".join(other_args3)}')
    log.info('Exit.')
    return


if __name__ == '__main__':
    bf_main()
