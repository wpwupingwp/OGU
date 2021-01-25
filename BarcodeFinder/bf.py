#!/usr/bin/python3

import argparse
import logging
from os.path import basename, exists
from os.path import join as join_path

import numpy as np
from Bio import SeqIO

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

    return


if __name__ == '__main__':
    main()
