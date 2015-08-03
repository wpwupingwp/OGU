#!/usr/bin/python3

import datetime
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

def Parser():
    '''Base on annotations in genbank files to extract fragments from Chloroplast Genome Sequence.
    '''
    Taxon = int(Record.features[0].qualifiers['db_xref'][0][6:])
    Organism = Record.annotations['organism']
    Accession = Record.annotations['accessions'][0]
    Gene = []
    All = []
    Type = 'whole'
    Start = 1
    End = len(Record)
    Sequence = str(Record.seq)
    Name = Organism
    Strand = 1
    rec = [Taxon, Organism, Accession, Name, Type, Start, End, Strand, Sequence, Date]
    All.append(rec)
    for i in Record.features:
        if i.type == 'gene' and 'gene' in i.qualifiers:
            if i.location_operator!= 'join':
                Type = 'gene'
                Start = int(i.location.start)
                End = int(i.location.end)
                Sequence = str(Record.seq[Start:End])
                Name = str(i.qualifiers['gene'][0])
                Strand = str(i.location.strand)
                rec = [Taxon, Organism, Accession, Name, Type, Start, End, Strand, Sequence, Date]
            elif i.location_operator == 'join':
                Type = 'gene'
                Start = int(i.sub_features[0].location.start)
                End = int(i.sub_features[0].location.end)
                Name = str(i.qualifiers['gene'][0])
                Strand = str(i.location.strand)
                Sequence = ''
                rec = [Taxon, Organism, Accession, Name, Type, Start, End, Strand, Sequence, Date]
                Gene.append(rec)
                Start = int(i.sub_features[1].location.start)
                End = int(i.sub_features[1].location.end)
                Sequence = ''.join([str(Record.seq[Start:End]), str(Record.seq[Start:End])])
                rec = [Taxon, Organism, Accession, Name, Type, Start, End, Strand, Sequence, Date]
        elif i.type == 'CDS' and 'gene' in i.qualifiers:
            Type = 'cds'
            Start = int(i.location.start)
            End = int(i.location.end)
            Sequence = str(Record.seq[Start:End])
            Name = str(i.qualifiers['gene'][0]).replace(' ', '_')
            Strand = str(i.location.strand)
            rec = [Taxon, Organism, Accession, Name, Type, Start, End, Strand, Sequence, Date]
        elif i.type == 'tRNA' and 'gene' in i.qualifiers:
            Type = 'tRNA'
            Start = int(i.location.start)
            End = int(i.location.end)
            Sequence = str(Record.seq[Start:End])
            if len(Sequence)>= 100:
                Sequence = ''
            Name = str(i.qualifiers['gene'][0]).replace(' ', '_')
            Strand = str(i.location.strand)
            rec = [Taxon, Organism, Accession, Name, Type, Start, End, Strand, Sequence, Date]
        elif i.type == 'rRNA':
            Type = 'rRNA'
            Start = int(i.location.start)
            End = int(i.location.end)
            Sequence = str(Record.seq[Start:End])
            Name = str(i.qualifiers['product'][0]).replace(' ', '_')
            Strand = str(i.location.strand)
            rec = [Taxon, Organism, Accession, Name, Type, Start, End, Strand, Sequence, Date]
        elif i.type == 'exon' and 'gene' in i.qualifiers :
            Type = 'exon'
            Start = int(i.location.start)
            End = int(i.location.end)
            Sequence = str(Record.seq[Start:End])
            if 'number' in i.qualifiers:
                Name = '_'.join([str(i.qualifiers['gene'][0]), 'exon', str(i.qualifiers['number'][0])])
            else:
                Name = '_'.join([str(i.qualifiers['gene'][0]), 'exon'])
            Strand = int(i.location.strand)
            rec = [Taxon, Organism, Accession, Name, Type, Start, End, Strand, Sequence, Date]
        elif i.type == 'intron' and 'gene' in i.qualifiers:
            Type = 'intron'
            Start = int(i.location.start)
            End = int(i.location.end)
            Sequence = str(Record.seq[Start:End])
            Strand = str(i.location.strand)
            if 'number' in i.qualifiers:
                Name = '_'.join([str(i.qualifiers['gene'][0]), 'intron', str(i.qualifiers['number'][0])])
            else:
                Name = '_'.join([str(i.qualifiers['gene'][0]), 'intron'])
            rec = [Taxon, Organism, Accession, Name, Type, Start, End, Strand, Sequence, Date]

        All.append(rec)
    Gene.sort(key = lambda x:x[5])

    for i in range(len(Gene)-1):
        Type = 'spacer'
        This = Gene[i]
        Next = Gene[i+1]
        Tail = This[6]+1
        Head = Next[5]-1
        Sequence = str(Record.seq[Tail:Head])
        Name = '_'.join(['-'.join([This[3], Next[3]]), 'Spacer'])
        Strand = 0
        rec = [Taxon, Organism, Accession, Name, Type, Start, End, Strand, Sequence, Date]
        All.append(rec)
    All.extend(Gene)

    SeqDB.extend(All)

def InitSeq():
    '''Init Sequence Database.
    '''
    con = sqlite3.connect('./test/DB')
    cur = con.cursor()
    cur.execute('create table if not exists main (Taxon int, Organism text, Accession text, Name text, Type text, Head int, Tail int,  Strand text, Sequence text, Date text, ID integer PRIMARY KEY);')
    for row in SeqDB:
        if row[8]!= '':
            cur.execute('insert into main (Taxon, Organism, Accession, Name, Type, Head, Tail, Strand, Sequence, Date) values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?);', (row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9]))
    con.commit()
    cur.close()
    con.close()
    print('Done.\n')
    
def SeqBatchQuery():
    con = sqlite3.connect('./test/DB')
    cur = con.cursor()
    listfile = input('list file name:\n')
    Results = list()
    with open(listfile, 'r') as In:
        Organisms = In.read().split(sep='\n')
    cur.execute('create table if not exists tasklist (Name text);')
    for Organism in Organisms:
        cur.execute('insert into tasklist (Name) values (?);', (Organism, ))
    cur.execute('select Taxon, Organism, Name, Type, Strand, Sequence, Head from main where Organism in (select Name from tasklist) order by Head', (Organism))
    Result = cur.fetchall()
    cur.execute('drop table tasklist;')
    cur.close()
    con.close()
    All = []
    for i in Result:
        Title = '|'.join([str(i[0]), i[1], i[2], i[3]])
        Filename = i[2]
        Sequence = MutableSeq(i[5])
        if i[4] == '-1':
            Sequence.seq = Sequence.reverse_complement()
        Record = [Title, Filename, Sequence]
        All.append(Record)
    for i in All:
        with open(''.join(['./out/', i[1], '.fasta']), 'a') as Fileout:
            Fileout.write('>%s\n%s\n'%(i[0], i[2]))
#rps12 may have larger than 50k fragments,  here to filter it
    rps12 = SeqIO.parse('./out/rps12.fasta', 'fasta')
    rps12short = list()
    for item in rps12:
        if len(item.seq)<4000:
            rps12short.append(item)
    SeqIO.write(rps12short, './out/rps12short.fasta', 'fasta')
    print('Done.\n')


def SeqQuery():
    '''Sequence query function,  to be continued.
    '''
    Querytype = input('1.Specific fragment\n2.Specific Organism\n3.Specific gene\n4.All\n5.All cds\n')
    organize = input('Organize output?(y/n)\n')
    if Querytype not in ['1', '2', '3', '4', '5']:
        raise ValueError('wrong input!\n')
    con = sqlite3.connect('./test/DB')
    cur = con.cursor()
    if Querytype == '1':
        Organism = input('Organism:\n')
        Gene = input('Gene:\n')
        Type = input('Fragment type(gene, cds, rRNA, tRNA, exon, intron, spacer):\n')
        cur.execute('select Taxon, Organism, Name, Type, Strand, Sequence from main where Name like ? and Type = ? and Organism=?', ('%'+Gene+'%', Type, Organism))
        Result = cur.fetchall()
    elif Querytype == '2':
        Organism = input('Organism:\n')
        Type = input('Fragment type(gene, cds, rRNA, tRNA, exon, intron, spacer, whole, fragments):\n')
        if Type == 'fragments':
            cur.execute('select Taxon, Organism, Name, Type, Strand, Sequence, Head from main where Organism = ?  order by Head', (Organism, ))
        else:
            cur.execute('select Taxon, Organism, Name, Type, Strand, Sequence, Head from main where Organism like ? and Type = ? order by Head', ('%'+Organism+'%', Type))
        Result = cur.fetchall()
    elif Querytype == '3':
        Gene = input('Gene:\n')
        Type = input('Fragment type(gene, cds, rRNA, tRNA, exon, intron, spacer):\n')
        cur.execute('select Taxon, Organism, Name, Type, Strand, Sequence from main where Name like ? and Type = ? order by Taxon', ('%'+Gene+'%', Type))
        Result = cur.fetchall()
    elif Querytype == '4':
        cur.execute('select Taxon, Organism, Name, Type, Strand, Sequence, Head from main order by Taxon')
        Result = cur.fetchall()
    elif Querytype == '5':
        cur.execute('select Taxon, Organism, Name, Type, Strand, Sequence, Head from main where type = "cds" order by Taxon')
        Result = cur.fetchall()

    All = []
    for i in Result:
        Title = '|'.join([str(i[0]), i[1], i[2], i[3]])
        Sequence = MutableSeq(i[5])
        gene = i[2]
        if i[4] == '-1':
            Sequence.seq = Sequence.reverse_complement()
        Record = [Title, gene, Sequence]
        All.append(Record)

    if organize ==  'y':
        makedirs('output')
        for i in All:
            file_name = ''.join([
                'output',
                '/',
                i[1].replace('/', ''),
                '.fasta'
            ])
            with open(file_name, 'a') as output_file:
                output_file.write('>%s\n%s\n' %(i[0], i[2]))
    else:
        Output = input('Enter output filename:\n')
        with  open('.'.join([Output, 'fasta']), 'w') as output_file:
            for i in All:
                output_file.write('>%s\n%s\n'%(i[0], i[2]))

    cur.close()
    con.close()
    print('Done.\n')

def UpdateSeqDBFromGenbank():
    '''Update Sequence database from Genbank,  need time to download.
    '''
    Down = urllib.request.urlopen('http://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=2759&opt=plastid').read().decode('utf-8')
    GenomeList = re.findall('((?<=nuccore/)[0-9]{9})', Down)
    userEmail = input('Input your email address for downloading data or use default(by press enter):\n')
    if userEmail is '\n':
        Entrez.email = 'wpwupingwp@outlook.com'
    else:
        Entrez.email = userEmail
    #need email address certify
    handle = Entrez.read(Entrez.epost(db='nuccore', id=', '.join(GenomeList)))
    W = handle['WebEnv']
    K = handle['QueryKey']
    GenomeContent = Entrez.efetch(db='nuccore', webenv=W, query_key=K, rettype='gb', retmode='text')
    Output = open('genbank', 'w')
    Output.write(GenomeContent.read())
    Output.close()
    UpdateSeqFromFile('genbank')
    
def UpdateSeqFromFile(FileIn):
    '''Update Sequence database from private file.
    '''
    global Record
    global SeqDB
    SeqDB = []
    handle = open(FileIn, 'r')
    Records = SeqIO.parse(FileIn, 'genbank')
    for Record in Records:
        Parser()
    InitSeq()
    handle.close()

def InitTaxon():
    '''Init Taxon database from file. 
    to be continued(add download function
    ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.tar.gz
    '''
    if exists('./test/taxdmp.zip') == False:
        ftp = FTP('ftp.ncbi.nih.gov')
        ftp.login()
        ftp.cwd('pub/taxonomy')
        ftp.retrbinary('RETR taxdump.zip', open('taxdmp.zip', 'wb').write)
        ftp.quit
    with ZipFile('./test/taxdmp.zip', 'r') as dmpfile:
        dmpfile.extractall(path='./test/')
    Id = dict()
    Data = list()
    Name = dict()
    Specie = list()
    Son = dict()
    GreatSon = dict()
    Parent = dict()
    Rank = dict()
    global Taxon
    Taxon = list()
    with open('./test/names.dmp', 'r') as dmpfile:
        raw = dmpfile.read().split(sep='\n')
        raw.pop()
        for record in raw:
            add = record.replace('\t', '').split(sep='|')
            if add[0] not in Name or add[2] == 'scientific name':
                Name[add[0]] = add[1]
    with open('./test/nodes.dmp', 'r') as dmpfile:
        raw = dmpfile.read().split(sep='\n')
        raw.pop()
        for record in raw:
            add = record.replace('\t', '').split(sep='|')
#1696063|Sarcocystis corvusi||scientific name|   
            Id[add[0]] = add[1]
            Rank[add[0]] = add[3]
            if add[2] == 'species':
                Specie.append(add[0])
    for specie in Specie:
        record = [specie, ]
        while Id[specie]!= '1' :
            record.append(Id[specie])
            specie = Id[specie]
#        if '33090' in record:
#            record.pop()
#            record.pop()
            Data.append(record)
    for data in Data:
        for n in range(len(data)):
            if data[n] not in Parent:
                Parent[data[n]] = data[(n+1):]
            if n == 0:
                continue
            if data[n] not in Son:
                Son[data[n]] = {data[n-1], }
            else:
                Son[data[n]].add(data[n-1])
            if data[n] not in GreatSon:
                GreatSon[data[n]] = {data[0], }
            else:
                GreatSon[data[n]].add(data[0])
    for specie in Name.items():
        if specie[0] not in Son:
            Son[specie[0]] = set()
        if specie[0] not in Parent:
            Parent[specie[0]] = list()
        if specie[0] not in GreatSon:
            GreatSon[specie[0]] = set()
        record = [specie[0], Name[specie[0]], Rank[specie[0]], Son[specie[0]], Parent[specie[0]], GreatSon[specie[0]]]
        Taxon.append(record)

    con = sqlite3.connect('./test/DB')
    cur = con.cursor()
    cur.execute('create table if not exists taxon (Id text, Name text, Rank text, Son text, Parent text, GreatSon text);')
    for line in Taxon:
        Son = ' '.join(line[3])
        Parent = ' '.join(line[4])
        GreatSon = ' '.join(line[5])
        cur.execute('insert into taxon (Id, Name, Rank, Son, Parent, GreatSon) values (?, ?, ?, ?, ?, ?);', (line[0], line[1], line[2], Son, Parent, GreatSon))
    con.commit()
    cur.close()
    con.close()
    print('Done.\n')
    
def TaxonQueryAuto(Name):
    '''Taxon query for seqquery,  may be remove.
    '''
    con = sqlite3.connect('./test/DB')
    cur = con.cursor()
    cur.execute('select Parent from taxon where  Name like ?;', ('%'+Name+'%', ))
    Result = cur.fetchall()
    Rank = dict()
    for record in Record:
        Rank[result[0]] = result[1]
    '''to be contiuned'''

def TaxonQueryNoAuto():
    '''Interactive query taxon database.
    '''
    while True:
        Querytype = input('1.by id\n2.by name\n')
        if Querytype not in ['1', '2']:
            return
        con = sqlite3.connect('./test/DB')
        cur = con.cursor()
        if Querytype == '1':
            Id = input('taxon id:\n')
            cur.execute('select * from taxon where Id = ?;', (Id, ))
            Result = cur.fetchall()
        elif Querytype == '2':
            Name = input('scientific name:\n')
#            cur.execute('select * from taxon where Name like ?;', ('%'+Name+'%', ))
            cur.execute('select * from taxon where Name = ?;', (Name, ))
            Result = cur.fetchall()
        cur.execute('select Id, Name from taxon;')
        Result2 = cur.fetchall()
        cur.close()
        con.close()
        Namedict = {'':''}
        for item in Result2:
            Namedict[item[0]] = item[1]
        for i in Result:
            Id = i[0]
            Name = i[1]
            Rank = i[2]
            Son = i[3].split(sep=' ')
            Sonname = list()
            for item in Son:
                Sonname.append(Namedict[item])
            Parent = i[4].split(sep=' ')
            Parentname = list()
            for item2 in Parent:
                Parentname.append(Namedict[item2])
            GreatSon = i[5].split(sep=' ')
            GreatSonname = list()
            for item3 in GreatSon:
                GreatSonname.append(Namedict[item3])
            handle = open('out.txt', 'a', encoding='utf-8')
            handle.write(''.join(['id      : ', Id, '\n']))
            handle.write(''.join(['name    : ', Name, '\n']))
            handle.write(''.join(['rank    : ', Rank, '\n']))
            handle.write(''.join(['parent  : ', '->'.join(Parentname), '\n']))
            handle.write(''.join(['son     : ', ', '.join(Sonname), '\n']))
            handle.write(''.join(['greatson: ', ', '.join(GreatSonname), '\n\n']))

#__main()__

'''main function,  entrance of the program.'''

Option = input('Select:\n1.Update database from GenBank\n2.Add pvirate data\n3.Query\n4.Init Taxon\n5.Query Taxon\n6.Batch Query\n')
Date = str(datetime.datetime.now())[:19].replace(' ', '-')
if Option == '1':
    UpdateSeqDBFromGenbank()
elif Option == '2':
    FileIn = input('Genbank format filename:\n')
    UpdateSeqFromFile(FileIn)
elif Option == '3':
    SeqQuery()
elif Option == '4':
    InitTaxon()
elif Option == '5':
    TaxonQueryNoAuto()
elif Option == '6':
    SeqBatchQuery()
else:
    print('Input error!\n')
