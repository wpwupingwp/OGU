#!/usr/bin/python3

import setuptools
from urllib.request import urlopen
from pathlib import Path
from zipfile import ZipFile


def init_lineage():
    """
    Only called by setup.py
    """
    data_folder = Path('data')
    url = 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip'
    with open('taxdmp.zip', 'wb') as out:
        out.write(urlopen(url).read())
    with ZipFile('taxdmp.zip', 'r') as dumpfile:
        dumpfile.extract('names.dmp', path='.')
        dumpfile.extract('nodes.dmp', path='.')
    id_name = {}
    id_rank = {}
    superkingdoms = []
    kingdoms = []
    phyla = []
    classes = []
    animal_orders = []
    with open('names.dmp', 'r') as names:
        for _ in names:
            line = _.split(sep='|')
            taxon_id = line[0].strip()
            name = line[1].strip()
            name_type = line[3].strip()
            if name_type == 'scientific name':
                id_name[taxon_id] = name
    with open('nodes.dmp', 'r') as nodes:
        for _ in nodes:
            line = _.split(sep='|')
            taxon_id = line[0].strip()
            rank = line[2].strip()
            id_rank[taxon_id] = rank
    for taxon_id in id_name:
        rank_name = id_rank[taxon_id]
        taxon_name = id_name[taxon_id]
        if rank_name == 'superkingdom':
            superkingdoms.append(taxon_name)
        elif rank_name == 'kingdom':
            kingdoms.append(taxon_name)
        elif rank_name == 'phylum':
            phyla.append(taxon_name)
        elif rank_name == 'class':
            classes.append(taxon_name)
        elif rank_name == 'order':
            if not taxon_name.startswith('ales'):
                animal_orders.append(taxon_name)

    with open(data_folder / 'superkingdoms.csv', 'w') as out:
        out.write(','.join(superkingdoms))
    with open(data_folder / 'kingdoms.csv', 'w') as out:
        out.write(','.join(kingdoms))
    with open(data_folder / 'phyla.csv', 'w') as out:
        out.write(','.join(phyla))
    with open(data_folder / 'classes.csv', 'w') as out:
        out.write(','.join(classes))
    with open(data_folder / 'animal_orders.csv', 'w') as out:
        out.write(','.join(animal_orders))
    return


init_lineage()

with open('README.md', 'r') as _:
    long_description = _.read()

with open('requirements.txt', 'r') as _:
    requires = [i.strip() for i in _.readlines()]

setuptools.setup(
    author='Ping Wu',
    author_email='wpwupingwp@outlook.com',
    description='All-in-one solution for discovering novel DNA barcode',
    install_requires=requires,
    include_package_data=True,
    # package_data={'BarcodeFinder': ['data/animal_orders.csv',
    #                                 'data/classes.csv', 'data/kingdoms.csv',
    #                                 'data/phyla.csv',
    #                                 'data/superkingdoms.csv']},
    license='GNU AGPL v3',
    long_description=long_description,
    long_description_content_type='text/markdown',
    name='BarcodeFinder',
    packages=setuptools.find_packages(),
    url='https://github.com/wpwupingwp/BarcodeFinder',
    version='0.9.44',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Affero General Public License v3',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
)
