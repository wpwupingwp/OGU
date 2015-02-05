#!/usr/bin/python3

import sys
import urllib.request
import json
import os.path
from Bio import SearchIO
from Bio.Blast import NCBIWWW

import warnings
from Bio import BiopythonExperimentalWarning
warnings.simplefilter('ignore',BiopythonExperimentalWarning)
#Try to ignore this warning, but doesn't work


def Query():
    with open(sys.argv[1],'r') as In:
        records=In.read()
    blast_results=NCBIWWW.qblast(
        program='blastn',
        database='nr',
        sequence=records,
        alignments=0,
        hitlist_size=3,     #Only use first three hits to identify the organism name of the query sequence
        format_type='XML'
    )
    with open(tmpfile_name,'w') as tmpfile:
        tmpfile.write(blast_results.read())

def Parse():
    parse_blast_results=list(SearchIO.parse(tmpfile_name,'blast-xml'))
    for result in parse_blast_results:
        query_id=result.description
        add=list()
        for hit in result:
            hit_id=hit.description
            hit_score=hit[0].evalue     #Only use first HSP
            name=hit_id.split(sep=' ')[:2]      #The first three words of description is usually enough to know the organism name
            name=' '.join(name)
            dictionary[name]=None
            add.append([query_id,hit_id,str(hit_score),name])
        out.extend(add)

def Translate():
    api_key='1630771459' #1000 times per hour
    for words in dictionary.keys():
        if dictionary[words] is not None:
            continue
        youdao_results=urllib.request.urlopen(''.join([
            'http://fanyi.youdao.com/openapi.do?keyfrom=Blastit&key=',
            api_key,
            '&type=data&doctype=json&version=1.1&q=',
            words
        ])).read()
        parse_translate_results=json.loads(youdao_results)
        translation_results=parse_translate_results['translation']
        dictionary[words]=translation_results[0]

def Output():
    n=0
    for item in out:
        if n%3 is 0:
            print('\nQuery sequence id:\t',item[0])
        print(
            '\t','Description:',item[1],'\n',
            '\t','Evalue:',item[2],'\n',
            '\t','Possible name:',item[3],'\n',
            '\t','Chinese:', dictionary[item[3]],
            '\n'
              )
        n=n+1

tmpfile_name='.'.join([sys.argv[1],'tmp'])
out=list()
dictionary=dict()
if os.path.exists(tmpfile_name):
    print('Found old results.')
    Parse()
    Output()
else:
    Query()
    Parse()
    Output()
