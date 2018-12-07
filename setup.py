#!/usr/bin/python3

import setuptools

with open('README.rst', 'r') as _:
    long_description = _.read()

requires = ['biopython>=1.72',
            'certifi>=2018.11.29',
            'matplotlib>=3.0.0',
            'numpy>=1.15.2',
            'pip>=18',
            'primer3-py>=0.5.7',
            'wheel>=0.32.3']

setuptools.setup(
    author='Ping Wu',
    author_email='wpwupingwp@outlook.com',
    description='All-in-one solution for discovering novel DNA barcode',
    install_requires=requires,
    license='GNU AGPL v3',
    long_description=long_description,
    long_description_content_type='text/markdown',
    name='BarcodeFinder',
    packages=setuptools.find_packages(),
    py_modules=['BarcodeFinder'],
    url='https://github.com/wpwupingwp/BarcodeFinder',
    version='0.9.15',
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
