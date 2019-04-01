#!/usr/bin/python3

import setuptools

with open('README.rst', 'r') as _:
    long_description = _.read()

with open('requirements.txt', 'r') as _:
    requires = [i.strip() for i in _.readlines()]

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
    version='0.9.33',
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
