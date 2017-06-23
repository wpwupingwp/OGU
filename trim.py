from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from timeit import default_timer as timer
import argparse


def main():
    start = timer()
    arg = argparse.ArgumentParser()
    arg.add_argument('input', help='input fasta file')
    arg.add_argument('location', help='fragment to cut, "x-y"')
    arg.add_argument('-n', action='store_true',
                     help='use N to replace where to cut')
    arg = arg.parse_args()
    start, end = arg.location.split('-')
    start = int(start) - 1
    end = int(end)
    with open('{}.new.fasta'.format(arg.input), 'w') as output:
        for record in SeqIO.parse(arg.input, 'fasta'):
            if arg.n:
                new = record[0:start] + SeqRecord(
                    'N'*(end-start), id=record.id,
                    description=record.description) + record[end:]
            else:
                new = record[0:start] + record[end:]
            SeqIO.write(new, output, 'fasta')
    end = timer()
    print('Cost {:.3f} seconds.'.format(end-start))


if __name__ == '__main__':
    main()
