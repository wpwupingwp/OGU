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
example, "\*.fasta" (include qutation mark) means all fasta files in the
folder.
### Sequence ID
BarcodeFinder use uniform sequence ID in all fasta files it generated.
```
SeqName|Order|Family|Genus|Species|Accession|SpecimenID
# example
rbcL|Poales|Poaceae|Oryza|longistaminata|MF998442|TAN:GB60B-2014
```
The order of the seven fields is fixed. Each field was seperated by the
vertical bar ("|"). The space character (" ") was disallowed and it was
replaced by underscore ("\_").
Because of technical issue, the order rank may be empty for animals.
* SeqName

    It means the name of sequences. Usually it its the gene name. For spacer,
    it is "geneA_geneB" that use underscore ("\_") to connect two gene's name.

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

    The list of primer pairs. CSV (comma-separated values text) format. The _b_
    is the name of the fragment (usually gene or spacer).

    Its title:
    ```
    Score,Sequences,AvgProductLength,StdEV,MinProductLength,MaxProductLength,Coverage,Resolution,TreeValue,AvgTerminalBranchLen,Entropy,LeftSeq,LeftTm,LeftAvgBitscore,LeftAvgMismatch,RightSeq,RightTm,RightAvgBitscore,RightAvgMismatch,DeltaTm,AlnStart,AlnEnd,AvgSeqStart,AvgSeqEnd
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

    * AvgTerminalBranchLen

        The average of sum of terminal branch's length.
    * Entropy

        <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;E_{H}&space;=&space;\frac{-&space;\sum_{i=1}^{k}{p_{i}&space;\log(p_{i})}}{\log(k)}" title="E_{H} = \frac{- \sum_{i=1}^{k}{p_{i} \log(p_{i})}}{\log(k)}" />

        The Shannon equitability index of the sub-alignment. The value is
        between 0 and 1.
    * LeftSeq

        Sequence of forward primer. The direction is 5' to 3'.
    * LeftTm

        The melting temperature of forward primer. The unit is Celcius (°C).
    * LeftAvgBitscore

        The average of raw bitscore of forward primer which is calculated by
        BLAST.
    * LeftAvgMismatch

        The average number of mismatch bases of forward primer which is
        counted by BLAST.
    * RightSeq

        Sequence of reverse primer. The direction is 5' to 3'.
    * RightTm

        The melting temperature of reverse primer. The unit is Celcius (°C).
    * RightAvgBitscore

        The average of raw bitscore of reverse primer which is calculated by
        BLAST.
    * RightAvgMismatch

        The average number of mismatch bases of reverse primer which is
        counted by BLAST.
    * DeltaTm

        The difference of melting temperature of forward and reverse primer.
        High DeltaTm may result in failure of PCR amplified.
    * AlnStart

        The location of beginning of forward primer (5', leftmost of primer
        pairs) in the whole alignment.
    * AlnEnd

        The location of end of reverse primer (5', rightmost of primer pairs)
        in the whole alignment.
    * AvgSeqStart

        The average beginning of forward primer in the original sequences.
        *ONLY USED FOR DEBUG*.
    * AvgSeqEnd

        The average end of forward primer in the original sequences.
        *ONLY USED FOR DEBUG*.

    The primer pairs were sorted by *Score*. Since the score may not fully
    satisfying user's specific consideration, it is suggested to choose primer
    pairs manually if the first primer pair failed on PCR experiment.
* _b_.fasta.uniq.fastq

    The fastq format file of primer's sequence. It contains two sequences and
    the direction is 5' to 3'. The first is the forward primer and the second
    is the reverse primer. The quality of each base equal to its proportion of
    the column in the alignment. Note that it may contains amibiguous base if
    you did not disable it.
* _b_.fasta.uniq.pdf

    The PDF format of the figure of sliding-window scan result of the
    alignment.
* _b_.fasta.uniq.png

    The PNG format of the figure of sliding-window scan result of the
    alignment.
* _b_.fasta.uniq.variance.tsv

    The CSV format of sliding-window scan result. The *"Index"* means the
    location of the base in the alignment. Note that the value DO NOT means
    the variance of the column of the base but the fragment started from this
    column.

* _b_.expand.fasta.uniq.\*

    By default, BarcodeFinder expands the alignment to its upstream/downstream
    to try to find suitable primer pairs for amplifed whole length of the
    target fragment. The filename contains "expand" does not ensure the final
    primers located in upstream/downstream. It only means BarcodeFinder used
    "expanded" data. *Users can use "-expand 0" to skip it.*
* Log.txt

    The log file. Contains all information printed in screen.
* Options.json

    The JSON file stored all options user inputed. It may be used by set
    "-json filename" to reduce input for similar run of the program.
* Summary.csv

    The summary of all ".fasta.uniq.csv" files which only contains information
    about the variance of each fragment. The only new field, *GapRatio*, means
    the ratio of the gap ("-") in the alignment. Higher value means the
    sequences may be too variaty that the alignment is not reliable.
* _b_-groupby_name

    The folder contains *"undivided"* sequences and intermediate results.
    Actually it is "roughly divided" sequences. The original genbank file was
    firstly divided into different fasta files if the genbank record contained
    different content. Usualy one genbank record contains serveral annotated
    region (multiple genes, for example). If two records contains same series
    of annotation (same order), they were put into same fasta file. Each file
    contains the intact sequence in the related genbank record.

* _b_-groupby_gene

    The folder contains divided sequences and intermediate results. After
    divided step occured in _b_-groupby_name, BarcodeFinder then divided each
    *cluster* of genbank records to several fasta files that each file
    contains only one region (one locus, one gene one spacer or one
    misc_feature) of the annotation.

    For instance, a record in "rbcL.gb" file may also contains atpB gene's
    sequence. The "rbcL.fasta" file does not contains any upstream/downstream
    sequence (except for ".expand" files) and "atpB_rbcL.fasta" does not have
    even one base of atpB or rbcL gene but only the spacer (assume that the
    annotation is precise).

    User can skip this divided step by set "-no_divide" to use the whole
    sequence for analysis. Note that it DOES NOT skip the first step of the
    dividing.

    These two folders usually can be ignored. However, sometimes, user may
    utilize one of these intermediate result:
    * _b_.fasta

        The raw fasta file directly converted from genbank file.
    * _b_.fasta.uniq

        Redundant sequences were removed.
    * _b_.fasta.uniq.aln

        The alignment of the fasta file.
    * _b_.fasta.uniq.aln.candidate.fasta

        The candidate primers. It may contains thousands of records. Do not
        recommend to pay attention to it.
    * _b_.fasta.uniq.aln.candidate.fasta.fastq

        Still, the candidate primers. This time it has quality which equals to
        base's proportion in the column of the alignment.
    * _b_.fasta.uniq.consensus.fastq

        The fastq format of the consensus sequence of the alignment.  Note
        that it contains aligment gap ("-").  Although this may be the most
        useful file in the folder, it is NOT RECOMMENDED to directly use it as
        the consensus sequence of the alignment because the
        consensus-genrating algorithm were optimized for primer design that it
        may be different with the _"real"_ consensus.

### Options
#### Help
#### General options
optional arguments:
  -h, --help            show this help message and exit
  -aln ALN              aligned fasta files to analyze (default: None)
  -fasta FASTA          unaligned fasta format data to add (default: None)
  -gb GB                genbank files (default: None)
  -j JSON               configuration json file (default: None)
  -stop {1,2,3}         Stop after which step: 1. Download 2. Preprocess data
                        3. Analyze (default: 3)
  -out OUT              output directory (default: None)

Genbank:
  -email EMAIL          email address for querying Genbank (default: None)
  -gene GENE            gene name (default: None)
  -group {animals,plants,fungi,protists,bacteria,archaea,viruses}
                        Species kind (default: None)
  -min_len MIN_LEN      minium length (default: 100)
  -max_len MAX_LEN      maximum length (default: 10000)
  -molecular {DNA,RNA}  molecular type (default: None)
  -organelle {mitochondrion,plastid,chloroplast}
                        organelle type (default: None)
  -query QUERY          query text (default: None)
  -taxon TAXON          Taxonomy name (default: None)

Preprocess:
  -expand EXPAND        expand length of upstream/downstream (default: 200)
  -max_name_len MAX_NAME_LEN
                        maximum length of feature name (default: 50)
  -max_seq_len MAX_SEQ_LEN
                        maximum length of sequence (default: 20000)
  -no_divide            analyze whole sequence instead of divided fragment
                        (default: False)
  -rename               try to rename gene (default: False)
  -uniq {longest,random,first,no}
                        method to remove redundant sequences (default: first)

Evaluate:
  -f                    faster evaluate variance by omit tree_valueand
                        terminal branch length (default: False)
  -s STEP               step length for sliding-window scan (default: 50)

Primer:
  -a AMBIGUOUS_BASE_N   number of ambiguous bases (default: 4)
  -c COVERAGE           minium coverage of base and primer (default: 0.6)
  -m MISMATCH           maximum mismatch bases in primer (default: 4)
  -pmin MIN_PRIMER      minimum primer length (default: 18)
  -pmax MAX_PRIMER      maximum primer length (default: 24)
  -r RESOLUTION         minium resolution (default: 0.5)
  -t TOP_N              keep n primers for each high varient region (default:
                        1)
  -tmin MIN_PRODUCT     minimum product length(include primer) (default: 300)
  -tmax MAX_PRODUCT     maximum product length(include primer) (default: 500)

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
