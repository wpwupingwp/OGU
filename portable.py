#!/usr/bin/python3

from datetime import datetime
from ftplib import FTP
import re
import sqlite3
import urllib.request
import warnings

from Bio import BiopythonDeprecationWarning
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import MutableSeq
from os import makedirs
from os.path import exists
from zipfile import ZipFile

warnings.simplefilter('ignore', BiopythonDeprecationWarning)


def parser(raw_seq, date):
    """Base on annotations in genbank files to extract fragments from Chloroplast Genome Sequence.
    """
    taxon_id = int(raw_seq.features[0].qualifiers['db_xref'][0][6:])
    organism = raw_seq.annotations['organism']
    accession = raw_seq.annotations['accessions'][0]
    gene = []
    records = []
    frag_type = 'whole'
    begin = 1
    end = len(raw_seq)
    sequence = str(raw_seq.seq)
    name = organism
    strand = 1
    rec = [taxon_id, organism, accession, name, frag_type, begin, end, strand, sequence, date]
    records.append(rec)
    for i in raw_seq.features:
        if i.type == 'gene' and 'gene' in i.qualifiers:
            if i.location_operator != 'join':
                frag_type = 'gene'
                begin = int(i.location.start)
                end = int(i.location.end)
                sequence = str(raw_seq.seq[begin:end])
                name = str(i.qualifiers['gene'][0])
                strand = str(i.location.strand)
                rec = [taxon_id, organism, accession, name, frag_type, begin, end, strand, sequence, date]
            elif i.location_operator == 'join':
                frag_type = 'gene'
                begin = int(i.sub_features[0].location.start)
                end = int(i.sub_features[0].location.end)
                name = str(i.qualifiers['gene'][0])
                strand = str(i.location.strand)
                sequence = ''
                rec = [taxon_id, organism, accession, name, frag_type, begin, end, strand, sequence, date]
                gene.append(rec)
                begin = int(i.sub_features[1].location.start)
                end = int(i.sub_features[1].location.end)
                sequence = ''.join([str(raw_seq.seq[begin:end]), str(raw_seq.seq[begin:end])])
                rec = [taxon_id, organism, accession, name, frag_type, begin, end, strand, sequence, date]
        elif i.type == 'CDS' and 'gene' in i.qualifiers:
            frag_type = 'cds'
            begin = int(i.location.start)
            end = int(i.location.end)
            sequence = str(raw_seq.seq[begin:end])
            name = str(i.qualifiers['gene'][0]).replace(' ', '_')
            strand = str(i.location.strand)
            rec = [taxon_id, organism, accession, name, frag_type, begin, end, strand, sequence, date]
        elif i.type == 'tRNA' and 'gene' in i.qualifiers:
            frag_type = 'tRNA'
            begin = int(i.location.start)
            end = int(i.location.end)
            sequence = str(raw_seq.seq[begin:end])
            if len(sequence) >= 100:
                sequence = ''
            name = str(i.qualifiers['gene'][0]).replace(' ', '_')
            strand = str(i.location.strand)
            rec = [taxon_id, organism, accession, name, frag_type, begin, end, strand, sequence, date]
        elif i.type == 'rRNA':
            frag_type = 'rRNA'
            begin = int(i.location.start)
            end = int(i.location.end)
            sequence = str(raw_seq.seq[begin:end])
            name = str(i.qualifiers['product'][0]).replace(' ', '_')
            strand = str(i.location.strand)
            rec = [taxon_id, organism, accession, name, frag_type, begin, end, strand, sequence, date]
        elif i.type == 'exon' and 'gene' in i.qualifiers:
            frag_type = 'exon'
            begin = int(i.location.start)
            end = int(i.location.end)
            sequence = str(raw_seq.seq[begin:end])
            if 'number' in i.qualifiers:
                name = '{0}_exon_{1}'.format(i.qualifiers['gene'][0],
                                             i.qualifiers['number'][0])
            else:
                name = '{0}_exon'.format(i.qualifiers['gene'][0])
            strand = int(i.location.strand)
            rec = [taxon_id, organism, accession, name, frag_type, begin, end, strand, sequence, date]
        elif i.type == 'intron' and 'gene' in i.qualifiers:
            frag_type = 'intron'
            begin = int(i.location.start)
            end = int(i.location.end)
            sequence = str(raw_seq.seq[begin:end])
            strand = str(i.location.strand)
            if 'number' in i.qualifiers:
                name = '{0}_{1}_intron'.format(i.qualifiers['gene'][0],
                                               i.qualifiers['number'][0])
            else:
                name = '{0}_intron'.format(i.qualifiers['gene'][0])
            rec = [taxon_id, organism, accession, name, frag_type, begin, end, strand, sequence, date]

        records.append(rec)
    gene.sort(key=lambda x: x[5])

    for i in range(len(gene) - 1):
        frag_type = 'spacer'
        now = gene[i]
        then = gene[i + 1]
        tail = now[6] + 1
        head = then[5] - 1
        sequence = str(raw_seq.seq[tail:head])
        name = '{0}-{1}_spacer'.format(now[3], then[3])
        strand = 0
        rec = [taxon_id, organism, accession, name, frag_type, begin, end, strand, sequence, date]
        records.append(rec)
    records.extend(gene)

    database.extend(records)


def init_seq():
    """Init Sequence Database.
    """
    con = sqlite3.connect('./data/DB')
    cur = con.cursor()
    cur.execute(
        'CREATE TABLE IF NOT EXISTS main (Taxon INT, Organism TEXT, Accession TEXT, Name TEXT, Type TEXT, Head INT, Tail INT,  Strand TEXT, Sequence TEXT, Date TEXT, ID INTEGER PRIMARY KEY);')
    for row in database:
        if row[8] != '':
            cur.execute(
                'INSERT INTO main (Taxon, Organism, Accession, Name, Type, Head, Tail, Strand, Sequence, Date) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?);',
                (row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9]))
    con.commit()
    cur.close()
    con.close()
    print('Done.\n')


def seq_batch_query():
    con = sqlite3.connect('./data/DB')
    cur = con.cursor()
    list_file = input('list file name:\n')
    with open(list_file, 'r') as In:
        organism_list = In.read().split(sep='\n')
    cur.execute('CREATE TABLE IF NOT EXISTS tasklist (Name TEXT);')
    for organism in organism_list:
        cur.execute('INSERT INTO tasklist (Name) VALUES (?);', (organism,))
    cur.execute(
        'SELECT Taxon, Organism, Name, Type, Strand, Sequence, Head FROM main WHERE Organism IN (SELECT Name FROM tasklist) ORDER BY Head',
        (organism))
    result = cur.fetchall()
    cur.execute('DROP TABLE tasklist;')
    cur.close()
    con.close()
    query_result = []
    for i in result:
        title = '{0}|{1}|{2}|{3}'.format(i[0], i[1], i[2], i[3])
        filename = i[2]
        sequence = MutableSeq(i[5])
        if i[4] == '-1':
            sequence.seq = sequence.reverse_complement()
        record = [title, filename, sequence]
        query_result.append(record)
    for i in query_result:
        with open('./out/{0}.fasta'.format(i[1]), 'a') as output_file:
            output_file.write('>{0}\n{1}\n'.format(i[0], i[2]))
            # rps12 may have larger than 50k fragments,  here to filter it
    rps12 = SeqIO.parse('./out/rps12.fasta', 'fasta')
    rps12short = list()
    for item in rps12:
        if len(item.seq) < 4000:
            rps12short.append(item)
    SeqIO.write(rps12short, './out/rps12short.fasta', 'fasta')
    print('Done.\n')


def seq_query():
    """Sequence query function,  to be continued.
    """
    query_type = input(
        '1.Specific fragment\n'
        '2.Specific Organism\n'
        '3.Specific gene\n'
        '4.All\n'
        '5.All cds\n'
    )
    organize = input('Organize output?(y/n)\n')
    if query_type not in ['1', '2', '3', '4', '5']:
        raise ValueError('wrong input!\n')
    con = sqlite3.connect('./data/DB')
    cur = con.cursor()
    if query_type == '1':
        organism = input('Organism:\n')
        gene = input('Gene:\n')
        frag_type = input('Fragment type(gene, cds, rRNA, tRNA, exon, intron, spacer):\n')
        cur.execute(
            'SELECT Taxon, Organism, Name, Type, Strand, Sequence FROM main WHERE Name LIKE ? AND Type = ? AND Organism=?',
            ('%' + gene + '%', frag_type, organism))
        result = cur.fetchall()
    elif query_type == '2':
        organism = input('Organism:\n')
        frag_type = input('Fragment type(gene, cds, rRNA, tRNA, exon, intron, spacer, whole, fragments):\n')
        if frag_type == 'fragments':
            cur.execute(
                'SELECT Taxon, Organism, Name, Type, Strand, Sequence, Head FROM main WHERE Organism = ?  ORDER BY Head',
                (organism,))
        else:
            cur.execute(
                'SELECT Taxon, Organism, Name, Type, Strand, Sequence, Head FROM main WHERE Organism LIKE ? AND Type = ? ORDER BY Head',
                ('%' + organism + '%', frag_type))
        result = cur.fetchall()
    elif query_type == '3':
        gene = input('Gene:\n')
        frag_type = input('Fragment type(gene, cds, rRNA, tRNA, exon, intron, spacer):\n')
        cur.execute(
            'SELECT Taxon, Organism, Name, Type, Strand, Sequence FROM main WHERE Name LIKE ? AND Type = ? ORDER BY Taxon',
            ('%' + gene + '%', frag_type))
        result = cur.fetchall()
    elif query_type == '4':
        cur.execute('SELECT Taxon, Organism, Name, Type, Strand, Sequence, Head FROM main ORDER BY Taxon')
        result = cur.fetchall()
    elif query_type == '5':
        cur.execute(
            'SELECT Taxon, Organism, Name, Type, Strand, Sequence, Head FROM main WHERE type = "cds" ORDER BY Taxon')
        result = cur.fetchall()

    query_result = []
    for i in result:
        title = '{0}|{1}|{2}|{3}'.format(i[0], i[1], i[2], i[3])
        sequence = MutableSeq(i[5])
        gene = i[2]
        if i[4] == '-1':
            sequence.seq = sequence.reverse_complement()
        record = [title, gene, sequence]
        query_result.append(record)

    if organize == 'y':
        if not exists('output'):
            makedirs('output')
        for i in query_result:
            file_name = 'output/{0}.fasta'.format(i[1].replace('/', ''))
            with open(file_name, 'a') as output_file:
                output_file.write('>{0}\n{1}\n'.format(i[0], i[2]))
    else:
        output = input('Enter output filename:\n')
        with open('{0}.fasta'.format(output), 'w') as output_file:
            for i in query_result:
                output_file.write('>{0}\n{1}\n'.format(i[0], i[2]))

    cur.close()
    con.close()
    print('Done.\n')


def update_seq_db_from_genbank(date):
    """Update Sequence database from Genbank,  need time to download.
    """
    down = urllib.request.urlopen(
        'http://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=2759&opt=plastid').read().decode('utf-8')
    genome_list = re.findall('((?<=nuccore/)[0-9]{9})', down)
    user_email = input('Input your email address for downloading data or use default(by press enter):\n')
    if user_email is '\n':
        Entrez.email = 'wpwupingwp@outlook.com'
    else:
        Entrez.email = user_email
    # need email address certify
    handle = Entrez.read(Entrez.epost(db='nuccore', id=', '.join(genome_list)))
    w = handle['WebEnv']
    k = handle['QueryKey']
    genome_content = Entrez.efetch(db='nuccore', webenv=w, query_key=k, rettype='gb', retmode='text')
    output = open('genbank', 'w')
    output.write(genome_content.read())
    output.close()
    update_seq_from_file('genbank', date)


def update_seq_from_file(genbank_file, date):
    """Update Sequence database from private file.
    """
    global database
    database = []
    records = SeqIO.parse(genbank_file, 'genbank')
    for raw_seq in records:
        parser(raw_seq, date)
    init_seq()


def init_taxon():
    """Init Taxon database from file.
    to be continued(add download function
    ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.tar.gz
    """
    if not exists('./data/taxdmp.zip'):
        ftp = FTP('ftp.ncbi.nih.gov')
        ftp.login()
        ftp.cwd('pub/taxonomy')
        ftp.retrbinary('RETR taxdmp.zip', open('./data/taxdmp.zip', 'wb').write)
        ftp.quit
    with ZipFile('./data/taxdmp.zip', 'r') as dumpfile:
        dumpfile.extractall(path='./data/')
    taxon_id = dict()
    data = list()
    name = dict()
    specie = list()
    son = dict()
    greatson = dict()
    parent = dict()
    rank = dict()
    global taxon
    taxon = list()
    with open('./data/names.dmp', 'r') as dumpfile:
        raw = dumpfile.read().split(sep='\n')
        raw.pop()
        for record in raw:
            add = record.replace('\t', '').split(sep='|')
            if add[0] not in name or add[2] == 'scientific name':
                name[add[0]] = add[1]
    with open('./data/nodes.dmp', 'r') as dumpfile:
        raw = dumpfile.read().split(sep='\n')
        raw.pop()
        for record in raw:
            add = record.replace('\t', '').split(sep='|')
            # 1696063|Sarcocystis corvusi||scientific name|
            taxon_id[add[0]] = add[1]
            rank[add[0]] = add[3]
            if add[2] == 'species':
                specie.append(add[0])
    for specie in specie:
        record = [specie, ]
        while taxon_id[specie] != '1':
            record.append(taxon_id[specie])
            specie = taxon_id[specie]
            #        if '33090' in record:
            #            record.pop()
            #            record.pop()
            data.append(record)
    for data in data:
        for n in range(len(data)):
            if data[n] not in parent:
                parent[data[n]] = data[(n + 1):]
            if n == 0:
                continue
            if data[n] not in son:
                son[data[n]] = {data[n - 1], }
            else:
                son[data[n]].add(data[n - 1])
            if data[n] not in greatson:
                greatson[data[n]] = {data[0], }
            else:
                greatson[data[n]].add(data[0])
    for specie in name.items():
        if specie[0] not in son:
            son[specie[0]] = set()
        if specie[0] not in parent:
            parent[specie[0]] = list()
        if specie[0] not in greatson:
            greatson[specie[0]] = set()
        record = [specie[0], name[specie[0]], rank[specie[0]], son[specie[0]], parent[specie[0]], greatson[specie[0]]]
        taxon.append(record)

    con = sqlite3.connect('./data/DB')
    cur = con.cursor()
    cur.execute(
        'CREATE TABLE IF NOT EXISTS taxon (Id TEXT, Name TEXT, Rank TEXT, Son TEXT, Parent TEXT, GreatSon TEXT);')
    for line in taxon:
        son = ' '.join(line[3])
        parent = ' '.join(line[4])
        greatson = ' '.join(line[5])
        cur.execute('INSERT INTO taxon (Id, Name, Rank, Son, Parent, GreatSon) VALUES (?, ?, ?, ?, ?, ?);',
                    (line[0], line[1], line[2], son, parent, greatson))
    con.commit()
    cur.close()
    con.close()
    print('Done.\n')


def taxon_query_auto(name):
    """Taxon query for seqquery,  may be remove.
     """
    con = sqlite3.connect('./data/DB')
    cur = con.cursor()
    cur.execute('SELECT Parent FROM taxon WHERE  Name LIKE ?;', ('%' + name + '%',))
    # Result = cur.fetchall()
    # Rank = dict()
    # for record in Record:
    #   Rank[result[0]] = result[1]
    '''to be continue'''


def taxon_query_no_auto():
    """Interactive query taxon database.
    """
    while True:
        query_type = input(
            '1.by id\n'
            '2.by name\n'
        )
        if query_type not in ['1', '2']:
            return
        con = sqlite3.connect('./data/DB')
        cur = con.cursor()
        if query_type == '1':
            taxon_id = input('taxon id:\n')
            cur.execute('SELECT * FROM taxon WHERE Id = ?;', (taxon_id,))
            result = cur.fetchall()
        elif query_type == '2':
            name = input('scientific name:\n')
            #            cur.execute('select * from taxon where Name like ?;', ('%'+Name+'%', ))
            cur.execute('SELECT * FROM taxon WHERE Name = ?;', (name,))
            result = cur.fetchall()
        cur.execute('SELECT Id, Name FROM taxon;')
        result2 = cur.fetchall()
        cur.close()
        con.close()
        name_dict = {'': ''}
        for item in result2:
            name_dict[item[0]] = item[1]
        for i in result:
            taxon_id = i[0]
            name = i[1]
            rank = i[2]
            son = i[3].split(sep=' ')
            son_name = list()
            for item in son:
                son_name.append(name_dict[item])
            parent = i[4].split(sep=' ')
            parent_name = list()
            for item2 in parent:
                parent_name.append(name_dict[item2])
            greatson = i[5].split(sep=' ')
            greatson_name = list()
            for item3 in greatson:
                greatson_name.append(name_dict[item3])
            handle = open('out.txt', 'a', encoding='utf-8')
            handle.write('id       : {0}\n'.format(taxon_id))
            handle.write('name     : {0}\n'.format(name))
            handle.write('rank     : {0}\n'.format(rank))
            handle.write('parent   : {0}\n'.format('->'.join(parent_name)))
            handle.write('son      : {0}\n'.format(', '.join(son_name)))
            handle.write('greatson : {0}\n\n'.format(', '.join(greatson_name)))


# __main()__
def main():
    '''Main function,  entrance of the program.'''
    option = input(
        'Select:\n'
        '1.Update database from GenBank\n'
        '2.Add pvirate data\n'
        '3.Query\n'
        '4.Init Taxon\n'
        '5.Query Taxon\n'
        '6.Batch Query\n'
    )
    date = str(datetime.now())[:19].replace(' ', '-')
    if option == '1':
        update_seq_db_from_genbank(date)
    elif option == '2':
        genbank_file = input('Genbank format filename:\n')
        update_seq_from_file(genbank_file, date)
    elif option == '3':
        seq_query()
    elif option == '4':
        init_taxon()
    elif option == '5':
        taxon_query_no_auto()
    elif option == '6':
        seq_batch_query()
    else:
        raise ValueError('Input error!\n')
    if not exists('data'):
        makedirs('data')

if __name__ == "__main__":
    main()
