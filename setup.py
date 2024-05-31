#!/usr/bin/python3
import setuptools
from OGU.global_vars import name
from OGU.utils import init_lineage

init_lineage()

with open('README.md', 'r') as _:
    long_description = _.read()

with open('requirements.txt', 'r') as _:
    requires = [i.strip() for i in _.readlines()]

setuptools.setup(
    author='Ping Wu',
    author_email='wpwupingwp@outlook.com',
    description='Organelle Genome Utilities',
    install_requires=requires,
    include_package_data=True,
    package_data={name: ['data/animal_orders.csv', 'data/classes.csv',
                         'data/kingdoms.csv', 'data/phyla.csv',
                         'data/superkingdoms.csv', 'data/button1.png',
                         'data/button2.png', 'data/button3.png',
                         'data/button4.png',
                         'visualize/draw_circle_plastid.ipynb',
                         'visualize/draw_circle_mitochondria.ipynb']},
    license='GNU AGPL v3',
    long_description=long_description,
    long_description_content_type='text/markdown',
    name=name,
    packages=setuptools.find_packages(),
    python_requires='>=3.8',
    entry_points={'console_scripts': ['OGU = OGU:main', 'ogu = OGU:main']},
    url='https://github.com/wpwupingwp/BarcodeFinder',
    version='1.60',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Affero General Public License v3',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
)
