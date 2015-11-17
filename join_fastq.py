from Bio import SeqIO
import sys

left = SeqIO.parse(sys.argv[1], 'fastq')
right = SeqIO.parse(sys.argv[2], 'fastq')
handle = open('combine.fastq', 'w')
offset = 64
for l,r in zip(left,right):
    print(l, r)
    l_seq = str(l.seq)
    r_seq = str(r.seq.reverse_complement())
    l_qual = ''.join([chr(i+offset) for i in
              l.letter_annotations['phred_quality']])
    r_qual = ''.join([chr(i+offset) for i in
              r.letter_annotations['phred_quality']])
    sequence = 'NNNNNNNNNN'.join([l_seq, r_seq])
    quality = 'AAAAAAAAAA'.join([l_qual, r_qual[::-1]])
    name = l.description
    handle.write('@{0}\n{1}\n+\n{2}'.format(name,sequence,quality))

