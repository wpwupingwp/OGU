#!/usr/bin/python3
import setuptools
from BarcodeFinder.global_vars import name
from BarcodeFinder.utils import init_lineage


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
    package_data={name: ['data/animal_orders.csv', 'data/classes.csv',
                                    'data/kingdoms.csv', 'data/phyla.csv',
                                    'data/superkingdoms.csv']},
    license='GNU AGPL v3',
    long_description=long_description,
    long_description_content_type='text/markdown',
    name=name,
    packages=setuptools.find_packages(),
    python_requires='>=3.6',
    url='https://github.com/wpwupingwp/BarcodeFinder',
    version='0.10.1',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Affero General Public License v3',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
)
