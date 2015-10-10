from glob import glob
from Bio import SeqIO

def prepare(filelist, result):
    for fasta_file in filelist:
        raw = list(SeqIO.parse(fasta_file, 'fasta'))
        for record in raw:
            length = len(record.seq)
            if length < min_length or length % 3 != 0:
                print('{0} in {1} has wrong length.\n'.format(record, 
                                                              fasta_file))

def count_gc(filelist,result):
            else:
                new_record = record[3:-3]
                gc1 = 0
                gc2 = 0
                gc3 = 0
                for i in range(len(new_record),3):
                    if new_record[i] in 'GC':
                        gc1 += 1
                    if new_record[i+1] in 'GC':
                        gc2 += 1
                    if new_record[i+2] in 'GC':
                        gc3 += 1
                result.append(','.join([record, fasta_file, gc1, gc2, gc3]))
    return 

def main():
    path = input('Enter the folder containing fasta files:\n')
    min_length = 300
    result = list()
    result.append([name, seq, gc1, gc2, gc3])
    filelist = glob(path)
    prepare(filelist, result)
    count_gc(result)
    with open('result.csv', 'w') as handle:
                  for i in result:
                      handle.write(','.join(i))
    return 0

if __name__ == '__main__':
    main()

