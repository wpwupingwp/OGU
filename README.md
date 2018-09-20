# BarcodeFinder
Automatic discover novel DNA barcode with universal primers.
## Prerequisite
### Software
* [Python3](https://www.python.org/downloads/)(3.5 or above)
* [BLAST](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
* [IQTREE](http://www.iqtree.org/)
* [MAFFT](https://mafft.cbrc.jp/alignment/software/)
### Python module
* biopython
* matplotlib
* numpy
* primer3-py
### Internet
The data retrive function requires Internet connection. Please ensure you have
stable network and inexpensive net fee when downloading large size of data.
## Installation
Assume that you alreadly installed python3 (3.5 or above), firstly, 
install BarcodeFinder.
```
# as administator
pip3 install BarcodeFinder
# normal user
pip3 install BarcodeFinder --user
```
Secondly, you need to install dependent software and python modules.
You can use install.py to help you:
```
# Windows
cd {BarcodeFinder path}
python install.py
# Linux and macos
cd {BarcodeFinder path}
python3 install.py
```
Here "BarcodeFinder path
##!@#$!@#$!
This program may help you to install all of them. Note that for Windows user,
you need to manually run BLAST installer and make sure that you added
installation path into _PATH_ enviroment variable..
## Usage
### Windows
python BarcodeFinder.py 
### Quick examples
Download all rbcL sequences of plants:
```
# Windows
python BarcodeFinder.py -query rbcL -group plants -stop 1 -out rbcL
# Linux and macos
python3 BarcodeFinder.py -query rbcL -group plants -stop 1 -out rbcL
```
Download all ITS sequences of Rosa and do preprocess:
```
# Windows
python BarcodeFinder.py -query "internal transcribed spacer" -taxon Rosa -stop 2 -out Rosa_its
# Linux and macos
python3 BarcodeFinder.py -query "internal transcribed spacer" -taxon Rosa -stop 2 -out Rosa_its
```
Download all Rosaceae chloroplast genome sequences, plus your data  as fasta 
format, then do analyze:
```
# Windows
python BarcodeFinder.py -organelle chloroplast -taxon Rosaceae -out Poaceae_cpg -fasta my_data.fasta
# Linux and macos
python3 BarcodeFinder.py -organelle chloroplast -taxon Rosaceae -out Poaceae_cpg -fasta my_data.fasta
```
Download sequences of Zea mays, set length between 100 bp and 3000 bp, plus
your aligned data, then do analyze:
```
# Windows
python BarcodeFinder.py -taxon "Zea mays" -min_len 100 -max_len 3000 -out Zea_mays -aln my_data.aln
# Linux and macos
python3 BarcodeFinder.py -taxon "Zea mays" -min_len 100 -max_len 3000 -out Zea_mays -aln my_data.aln
```
Download all Oryza chloroplast genomes, divide them into fragments, and skip
analyze:
```
# Windows
python BarcodeFinder.py -taxon Oryza -organelle chloroplast -stop 2 -out Oryza_cp
# Linux and macos
python3 BarcodeFinder.py -taxon Oryza -organelle chloroplast -stop 2 -out Oryza_cp
```

