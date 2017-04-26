#!/usr/bin/python3

from Bio import SeqIO
from sys import argv


handle = open(argv[1]+'.fasta', 'w')

for record in SeqIO.parse(argv[1], 'gb'):
    """
    From Zhang guojin
    order statswith ales
    family startswith aceae except 8
    http://duocet.ibiodiversity.net/index.php?title=%E4%BA%92%E7%94%A8%E5%90%8D
    %E7%A7%B0&mobileaction=toggle_view_mobile"""
    # order|family|organims(genus|species)
    order_family = record.annotations['taxonomy']
    family_exception_raw = (
        'Umbelliferae,Palmae,Compositae,Cruciferae,Guttiferae,Leguminosae,'
        'Leguminosae,Papilionaceae,Labiatae,Gramineae')
    family_exception = family_exception_raw[0].split(',')
    for item in order_family:
        if item.endswith('ales'):
            order = item
            continue
        if item.endswith('aceae'):
            family = item
        elif item in family_exception:
            family = item
    organism = record.annotations['organism'].replace(' ', '_')
    genus, *species = organism.split('_')
    taxon = '{}|{}|{}|{}'.format(order, family, genus, '_'.join(species))

    accession = record.annotations['accessions'][0]
    try:
        specimen = record.features[0].qualifiers['specimen_voucher'
                                                 ][0].replace(' ', '_')
    except:
        specimen = ''
    for feature in record.features:
        if feature.type == 'gene' and 'gene' in feature.qualifiers:
            gene = feature.qualifiers['gene'][0].replace(' ', '_')
            sequence = record.seq[feature.location.start:feature.location.end]
            sequence = str(sequence)
            handle.write('>{}|{}|{}|{}\n{}\n'.format(gene, taxon, accession,
                                                     specimen, sequence))
