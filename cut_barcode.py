from Bio import SeqIO
from time import process_time
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('-b', '--barcode_length', type=int, default=10, help='original barcode length')
parser.add_argument('-c', '--cut', type=int, default=5, help='barcode cutted length')
parser.add_argument('-l', '--linker', default='GTAGACTGCGTACC', help='linker length')
parser.add_argument('-f', '--barcode_file', help='barcode_file')
parser.add_argument('input', help='input fastq file')
parser.add_argument('-o', '--output', help='output fastq file')
parser.add_argument('-t', '--barcode_type',default='5*2')
arg = parser.parse_args()
if arg.barcode_type != '5*2':
    raise ValueError('Do not support this type of barcode.\n')

check_length = arg.barcode_length - arg.cut + len(arg.linker) - 1
with open(arg.barcode_file, 'r') as barcode_file:
    for i in barcode_file.readline():
        if i.startswith('>sample'):
            break
    barcode = re.findall('(?<=,)[ATCG]+(?=,)', barcode_file.read())
barcode_dict = {i:None for i in barcode}
barcode_link = dict()
for i in barcode:
    barcode_link[''.join([i[0:(arg.cut - 1 )], arg.linker])] = None

handle = open(arg.output, 'a') 
raw = SeqIO.parse(arg.input, 'fastq')
for sequence in raw:
    if str(sequence.seq[:arg.barcode_length]) in barcode_dict:
        pass
    elif str(sequence.seq[arg.cut:check_length]) in barcode_link:
        pass
    else:
        continue
    SeqIO.write(sequence[arg.cut:], handle, 'fastq')

print('Finished within {:.3f} seconds.'.format(process_time()))
