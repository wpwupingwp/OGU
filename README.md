# BarcodeFinder
BarcodeFinder could automatically discover novel DNA barcodes with universal
primers. It does three things as listed below.
* Collect data.
It can automatically retrieve data from NCBI Genbank with restrictions that
user provided, such as gene name, taxonomy, sequence name and organelle.
Also, it can integrate user provided sequences or alignments.
* Preprocess data
Barcodefinder utilizes annotation information in data to divide sequence into
fragments (gene, spacer, misc_feature), because data collected from Genbank
may not be "uniform". For instance, you can find one gene's upstream and
downstream sequences in one record but only gene sequence in another record.
The situation becomes worse for intergenic spacers, that various annotation
style may cause endless trouble in following analysis.
Given that one gene or spacer for one species may be sequenced several times,
by default, BarcodeFinder removes redundant sequences to left only one record
for each species. This behavior can be changed as you wish. Then, _mafft_ was
called for alignment. Each sequence's direction were adjusted and all
sequences were reordered.
* Analyze
Firstly, BarcodeFinder iterately evaluate variance of each alignment by
calculating Pi, Shannon Index, observed resolution, tree value and terminal
branch length, etc. If the result is lower than given threshold, i.e., it does
not have efficient resolution, this alignment were skipped.
Next, a sliding-window scan will be performed for those alignments passed the
test. The high-variance region (variance "hotspot") were picked and its
upstream/downstream region were used to find primer.
In those conserved region for finding primers, consensus sequences were
generated and with the help of primer3, candidate primers were selected.
After BLAST validation, suitable primers were combined to form serveral primer
pairs. According to the limit of PCR product's length, only pairs with wanted
length were left. Note that gaps were removed to calculated real length
instead of alignment length. The resolution of the subalignment were
recalculated to remove false positive primer pairs.
Finally, primer pairs were reorderd by score to make it easy for user to find
"best" primer pairs they want.
## Prerequisite
### Software
* Python3 (3.5 or above)
* BLAST+
* IQTREE
* MAFFT
### Python module
* Biopython
* matplotlib
* numpy
* primer3-py
### Internet
The data retrieve function requires Internet connection. Please ensure you have
stable network and inexpensive net fee when downloading large size of data.
## Installation
Assume that you already installed [Python3](https://www.python.org/downloads/)
(3.5 or above). Firstly, install BarcodeFinder.
The easiest way is to download
[BarcodeFinder.py](https://github.com/wpwupingwp/BarcodeFinder) and put it
into wherever you want.
If you would like to use pip, then:
```
# as administator
pip3 install BarcodeFinder
# normal user
pip3 install BarcodeFinder --user
```
Secondly, you need to install dependent software and python modules.
Although BarcodeFinder has assistant function to automatically install
dependent software and modules if it cannot find them. However, it is
highly recommended to follow official installation procedure to make it
easy for management and give you a clean working directory.
For Linux user, if you have root privileges, just use your package manager:
```
# Ubuntu and Debian
sudo apt install mafft ncbi-blast+ iqtree
# Fedora
sudo dnf install mafft ncbi-blast+ iqtree
# Fedora 2
sudo yum install mafft ncbi-blast+ iqtree
# ArchLinux
sudo pacman -S mafft ncbi-blast+ iqtree
# FreeBSD
sudo pkg install mafft ncbi-blast+ iqtree
```
For MacOS user, assume you have root privileges, if you do not have
*brew*, install it:
```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```
Then:
```
brew install blast mafft brewsci/science/iqtree
```
If you are Windows user or you do not have root privileges, follow these
instructions:
1. BLAST+
* [Windows](https://www.ncbi.nlm.nih.gov/books/NBK52637/)
* [Linux and MacOS](https://www.ncbi.nlm.nih.gov/books/NBK52640/)
2. MAFFT
* [Windows](https://mafft.cbrc.jp/alignment/software/windows.html)
Choose "All-in-one version", download and unzip. Then follow the step in
BLAST+ installation manual to set _PATH_.
* [Linux](https://mafft.cbrc.jp/alignment/software/linux.html)
Choose "Portable package", download and unzip. Then follow the instruction of
BLAST+ to set _PATH_.
* [MacOS](https://mafft.cbrc.jp/alignment/software/macosx.html)
Choose "All-in-one version", download and unzip. Then follow the step in
BLAST+ installation manual to set _PATH_.
3. IQTREE
[Download](http://www.iqtree.org/#download) according to your OS.
Unzip and add the path of subfolder *bin* into _PATH_
## Usage
The basic usage looks like this:
```
# Windows
python BarcodeFinder.py [input] -[options] -out [out_folder]
# Linux and MacOS
python3 BarcodeFinder.py [input] -[options] -out [out_folder]
```
### Input
BarcodeFinder accepts:
1. Genbank query. You can use "-query" or combine with other filters.
2. Unaligned fasta files. Each file were considered to be one locus to
   evaluate variance.
3. Alignments (fasta format).
For _2_ and _3_, ambiguous bases were allowed in sequence.
### Options

*[data]* means input. It can be Genbank query, fasta file names, alignments or
combinations. If you want to use "\*" or "?" to represent a series of files,
make sure to use _"_ to quote it. Also, if 

```
# Windows
python Barcodefinder.py -h
# Linux and MacOS
python3 Barcodefinder.py -h
```

## Output
BarcodeFinder renames all sequences in this model:
```
gene|order|family|genus|species|accesion_id
```
Here _gene_ means fragment's name.
The last thing in this step is to
## Quick examples
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
Download all Rosaceae chloroplast genome sequences, plus your data as fasta
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

