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
pip3 install BarcodeFinder
```
If you do not have administrator privilege, then use this instead:
```
pip3 install BarcodeFinder --user
```
If you want to update from old version:
```
pip3 install -U BarcodeFinder
```
Secondly, you need to install BLAST, IQTREE and MAFFT if you have not yet.
You can install them manually.
For Microsoft Windows users,
