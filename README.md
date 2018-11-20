# BarcodeFinder
BarcodeFinder could automatically discover novel DNA barcodes with universal
primers. It does three things as listed below.
* Collect data

    It can automatically retrieve data from NCBI Genbank with restrictions
    that user provided, such as gene name, taxonomy, sequence name and
    organelle. Also, it can integrate user provided sequences or alignments.
* Preprocess data

    Barcodefinder utilizes annotation information in data to divide sequence
    into fragments (gene, spacer, misc_feature), because data collected from
    Genbank may not be "uniform". For instance, you can find one gene's
    upstream and downstream sequences in one record but only gene sequence in
    another record. The situation becomes worse for intergenic spacers, that
    various annotation style may cause endless trouble in following analysis.

    Given that one gene or spacer for each species may be sequenced several
    times, by default, BarcodeFinder removes redundant sequences to left only
    one record for each species. This behavior can be changed as you wish.
    Then, _mafft_ was called for alignment. Each sequence's direction were
    adjusted and all sequences were reordered.
* Analyze

    Firstly, BarcodeFinder evaluate variance of each alignment by calculating
    Pi, Shannon Index, observed resolution, tree value and terminal branch
    length, etc. If the result is lower than given threshold, i.e., it does
    not have efficient resolution, this alignment were skipped.

    Next, a sliding-window scan will be performed for those alignments passed
    the test. The high-variance region (variance "hotspot") were picked and
    its upstream/downstream region were used to find primer.

    Consensus sequence of those conserved region for finding primers were
    generated and with the help of primer3, candidate primers were selected.
    After BLAST validation, suitable primers were combined to form several
    primer pairs. According to the limit of PCR product's length, only pairs
    with wanted length were left. Note that gaps were removed to calculated
    real length instead of alignment length. The resolution of the
    sub-alignment were recalculated to remove false positive primer pairs.

    Finally, primer pairs were reordered by score to make it easy for user to
    find "best" primer pairs they want.
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
The data retrieve function requires Internet connection. Please ensure you
have stable network and reasonable Internet traffic charge for downloading
large size of data.
## Installation
Assume that you already installed [Python3](https://www.python.org/downloads/)
(3.5 or above). Firstly, install BarcodeFinder.
The easiest way is to use pip:
```
# As administator
pip3 install BarcodeFinder
# Normal user
pip3 install BarcodeFinder --user
```
Secondly, you need to install dependent software and python modules.
BarcodeFinder has assistant function to automatically install dependent
software if it cannot find them. However, it is *highly recommended* to follow
official installation procedure to make it easy for management and give you a
clean working directory.

For Linux user, if you have root privileges, just use your package manager:
```
# Ubuntu and Debian
sudo apt install mafft ncbi-blast+ iqtree
# Fedora (1)
sudo dnf install mafft ncbi-blast+ iqtree
# Fedora (2)
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
    Choose "Portable package", download and unzip. Then follow the instruction
    of BLAST+ to set _PATH_.
    * [MacOS](https://mafft.cbrc.jp/alignment/software/macosx.html)
    Choose "All-in-one version", download and unzip. Then follow the step in
    BLAST+ installation manual to set _PATH_.
3. IQTREE
[Download](http://www.iqtree.org/#download)
Download installer according to your OS. Unzip and add the path of subfolder
*bin* into _PATH_
## Usage
BarcodeFinder is a command line program. Once you open your command line
(Windows) or
terminal (Linux and MacOS)
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
4. Genbank format files.
Note that ambiguous bases were allowed in sequence. If you want to use "\*" or
"?" to represent a series of files, make sure to use _"_ to quote it. For
example, "\*.fasta"(include qutation mark) means all fasta files in the
folder.
### Sequence ID
BarcodeFinder use uniform sequence ID in all fasta files it generated.
```
SeqName|Order|Family|Genus|Species|Accession|SpecimenID
# example
rbcL|Poales|Poaceae|Oryza|longistaminata|MF998442|TAN:GB60B-2014
```
The order of the seven fields is fixed. Each field was seperated by the
vertical bar ("|"). The space character(" ") was disallowed and it was
replaced by underscore("\_").
Because of technical issue, the order rank may be empty for animals.
* SeqName

    It means the name of sequences. Usually it its the gene name. For spacer,
    it is "geneA_geneB" that use underscore("\_") to connect two gene's name.

    Note that the original gene name in genbank file were renamed to try to
    fix part of annotation error. You can use raw name by set "-no_rename"
    option.

    If a valid sequence name could not be found in annotation of genbank file,
    BarcodeFinder will use "Unknown" instead.
* Order

    The order name of the species. 
* Family

    The family name of the species.
* Genus

    The genus name of the species, i.e., the first part of the scientific
    name.
* Species

    The specific epithet of the species, i.e., the second part of the
    scientific name of the species. It may contains subspecies name.
* Accession

    The Genbank Accession number of the sequence. Do not contain version
    number.
* SpecimenID

    The ID of specimen of the sequence. May be empty.

### Output
All results will be put in the output folder. If you didn't set output path by
"-out", BarcodeFinder will create a folder named by current time, for example, 
"2018-11-19T16-41-59.330217".
* _a_.gb
    
    The raw genbank file. The _a_ comes from the keyword of query.
* _a_.gb.fasta

    The converted fasta file of the ".gb" file.
* _b_.fasta.uniq.csv

    The list of primer pairs.
    Its title:
    ```
    Score,Sequences,AvgProductLength,StdEV,MinProductLength,MaxProductLength,Coverage,Resolution,TreeValue,Entropy,LeftSeq,LeftTm,LeftAvgBitscore,LeftAvgMismatch,RightSeq,RightTm,RightAvgBitscore,RightAvgMismatch,DeltaTm,AlnStart,AlnEnd,AvgSeqStart,AvgSeqEnd
    ```

    * Score

        The score of this primer pair. Usually the higher, the better.
    * Sequences
        
        How many sequences were used to find this primer pair.

    * AvgProductLength

        The average length of amplified DNA fragment by this pair of primer.
        Integer.
    * StdEV

        The standard deviation of the AvgProductLength. Higher number means
        the primer may amplified different length of DNA fragment. Lower
        number or even zero means close length, i.e., much more conservative.
    * MinProductLength

        The minimum length of amplified fragment.
    * MaxProductLength

        The maximum length of amplified fragment. Note that all these four
        fields were calculated by given sequences that may not cover
        exception.
    * Coverage

        The coverage of this primer pair for sequences it used. Calculated by
        BLAST result. High coverage means it is much more "universal".
    * Resolution

        <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;R_{o}=\frac{n_{uniq}}{n_{total}}" title="R_{o}=\frac{n_{uniq}}{n_{total}}" />

        The *observed resolution* of the sub-alignment sliced by the primer
        pair, which is equal to number of uniq sequences divided by number of
        total sequences. The value is
        between 0 and 1.
    * TreeValue

        <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;R_{T}=\frac{n_{internal}}{n_{terminal}}" title="R_{T}=\frac{n_{internal}}{n_{terminal}}" />

        The *tree resolution* of the sub-alignment, which is equal to number
        of internal nodes of phylogenetic tree construted from the alignment
        divided by number of terminal nodes. The value is between 0 and 1.
    * Entropy

        <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;E_{H}&space;=&space;\frac{-&space;\sum_{i=1}^{k}{p_{i}&space;\log(p_{i})}}{\log(k)}" title="E_{H} = \frac{- \sum_{i=1}^{k}{p_{i} \log(p_{i})}}{\log(k)}" />

        The Shannon equitability index of the sub-alignment. The value is
        between 0 and 1.
    * 


* _b_.fasta.uniq.resolution.tsv
* _b_.fasta.uniq.csv
* _b_.fasta.uniq.fastq
* _b_.fasta.uniq.pdf
* _b_.fasta.uniq.png
* _b_.expand.fasta.uniq.csv


# to be continue
BarcodeFinder renames all sequences in this model:
```
gene|order|family|genus|species|accesion_id
```
Here *gene* means fragment's name.
The last thing in this step is to
### Options
# to be continue
```
# Windows
python Barcodefinder.py -h
# Linux and MacOS
python3 Barcodefinder.py -h
```
## Quick examples
1. Download all _rbcL_ sequences of plants and do pre-process:
```
# Windows
python BarcodeFinder.py -query rbcL -group plants -stop 1 -out rbcL
# Linux and macos
python3 BarcodeFinder.py -query rbcL -group plants -stop 1 -out rbcL
```
2. Download all ITS sequences of _Rosa_ and do pre-process:
```
# Windows
python BarcodeFinder.py -query "internal transcribed spacer" -taxon Rosa -stop 2 -out Rosa_its
# Linux and macos
python3 BarcodeFinder.py -query "internal transcribed spacer" -taxon Rosa -stop 2 -out Rosa_its
```
3. Download all Rosaceae chloroplast genome sequences, plus your own data.
   Then do analyze:
```
# Windows
python BarcodeFinder.py -organelle chloroplast -taxon Rosaceae -out Poaceae_cpg -fasta my_data.fasta
# Linux and macos
python3 BarcodeFinder.py -organelle chloroplast -taxon Rosaceae -out Poaceae_cpg -fasta my_data.fasta
```
4. Download sequences of _Zea mays_, set length between 100 bp and 3000 bp,
   plus your aligned data, then do analyze:
```
# Windows
python BarcodeFinder.py -taxon "Zea mays" -min_len 100 -max_len 3000 -out Zea_mays -aln my_data.aln
# Linux and macos
python3 BarcodeFinder.py -taxon "Zea mays" -min_len 100 -max_len 3000 -out Zea_mays -aln my_data.aln
```
5. Download all _Oryza_ chloroplast genomes, divide them into fragments, and
   skip analyze:
```
# Windows
python BarcodeFinder.py -taxon Oryza -organelle chloroplast -stop 2 -out Oryza_cp
# Linux and macos
python3 BarcodeFinder.py -taxon Oryza -organelle chloroplast -stop 2 -out Oryza_cp
```
