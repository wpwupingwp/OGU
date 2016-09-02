#!/usr/bin/python3

import argparse
from Bio import SeqIO
from glob import glob
from time import process_time


def screen_length(fasta_files, cut_off):
    for fasta in fasta_files:
        keep = 0
        drop = 0
        count = list()
        handle = open(fasta, 'r')
        handle_out = open(fasta.replace('.fasta', '.fas'), 'w')
        f = SeqIO.parse(handle, 'fasta')
        for record in f:
            count.append(len(record))
        average = sum(count) / len(count)
        handle.seek(0)
        f = SeqIO.parse(handle, 'fasta')
        for record in f:
            print(len(record))
            if len(record) > average*cut_off:
                keep += SeqIO.write(record, handle_out, 'fasta')
            else:
                drop += 1
        print('{0} sequences longer than {1} were kept, dropped {2} short sequences.'.format(keep, average*cut_off, drop))
        handle.close()


def main():
    """This program will try to find out barcode to devide different species
    while ignore distinction among subspecies level."""
    parser = argparse.ArgumentParser( description=main.__doc__)
    parser.add_argument('--path', default='.', 
                        help='target path, default is present directory')
    parser.add_argument('--cut_off', default=0.8, type=int, help='cut off for screening length, those length shorter than average_length*cut_off will be drop off')
    parser.print_help()
    arg = parser.parse_args() 
    fasta_files = glob('*.fasta')
    screen_length(fasta_files, arg.cut_off)
    print('Cost {:.3f}s.\n'.format(process_time()))


if __name__ == '__main__':
    main()