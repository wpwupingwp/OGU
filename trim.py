from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from sys import argv

start, end = argv[2].split('-')
start = int(start) - 1
end = int(end)
with open('{}.new.fasta'.format(argv[1]), 'w') as output:
    for record in SeqIO.parse(argv[1], 'fasta'):
        # new = record[0:start] + record[end:]
        new = record[0:start] + SeqRecord('N'*(end-start), id=record.id,
                                          description = record.description) + record[end:]
        SeqIO.write(new, output, 'fasta')
