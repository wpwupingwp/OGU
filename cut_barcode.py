from Bio import SeqIO
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('-b', '--barcode_length', type=int, default=10, help='original barcode length')
parser.add_argument('-c', '--cut', type=int, default=5, help='barcode cutted length')
parser.add_argument('-l', '--linker', type=int, default=10, help='linker length')
parser.add_argument('f', 'barcode_file')
parser.add_argument('i', 'input', help='input fastq file')
parser.add_argument('o', 'output', help='output fastq file')
parser.add_argument('-t', '--barcode_type',default='5*2')
arg = parser.parse_args()
#print(arg)
if arg.t != '5*2':
    raise ValueError('Do not support this type of barcode.\n')

check_length = len(arg.b-arg.c+arg.l)-1
with open(arg.f, 'r') as barcode_file:
    barcode = re.findall('[ATCG]{0}(?>=,)'.format(arg.b), barcode_file.read())
barcode_dict = {i:None for i in barcode}
barcode_link = dict()
for i in barcode:
    barcode_link[''.join([i[0:(arg.cut - 1 )], arg.link])] = None

handle = open(arg.o, 'a') 
raw = SeqIO.parse(arg.i, 'fastq')
for sequence in raw:
    if str(sequence.seq[:arg.b]) in barcode_dict:
        pass
    elif str(sequence.seq[arg.cut:check_length]) in barcode_link:
        pass
    else:
        continue
    SeqIO.write(sequence[arg.cut:], handle, 'fastq')
