#!/usr/bin/python3

import argparse
import os
from Bio import SeqIO
from glob import glob
from time import process_time
from multiprocessing import cpu_count
from subprocess import run


def screen_merge(fasta_files, path, cut_off, sample):
    merge = list()
    output = list()
    for fasta in fasta_files:
        count = list()
        handle = open(fasta, 'r')
        file_name = fasta.replace('.fasta', '.fas')
        handle_out = open(file_name, 'w')
        output.append(file_name)
        f = SeqIO.parse(handle, 'fasta')
        for record in f:
            count.append(len(record))
        average = sum(count) / len(count)
        handle.seek(0)
        f = SeqIO.parse(handle, 'fasta')
        for n, record in zip(range(sample), f):
            if len(record) > average*cut_off:
                merge.append(record)
                SeqIO.write(record, handle_out, 'fasta')
        handle.close()
        handle_out.close()
    merge_file = path + '/merge.fas'
    SeqIO.write(merge, merge_file, 'fasta')
    output.append(merge_file)
    return output


def mafft(fas_file):
    for fas in fas_file:
        run(['mafft', '--reorder --thread {0} {1} > {2}'.format(cpu_count(), fas, fas.replace('.fas', '.aln'))], shell=True)


def main():
    """This program will try to find out barcode to devide different species
    while ignore distinction among subspecies level."""
    parser = argparse.ArgumentParser( description=main.__doc__)
    parser.add_argument('--path', default='.', 
                        help='target path, default is present directory')
    parser.add_argument('--cut_off', default=0.8, type=float, help='cut off for screening length, those length shorter than average_length*cut_off will be drop off')
    parser.add_argument('--sample',default=5,type=int,help='sample selected from each group')
    parser.print_help()
    arg = parser.parse_args() 
    fasta_files = glob(arg.path+'/*.fasta')
    #merge = screen_merge(fasta_files, arg.path, arg.cut_off, arg.sample)
    merge = screen_merge(fasta_files, **vars(arg))
    mafft(merge)
    print('Cost {:.3f}s.\n'.format(process_time()))


if __name__ == '__main__':
    main()