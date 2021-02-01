[![Build Status](https://travis-ci.org/wpwupingwp/novowrap.svg?branch=master)](https://travis-ci.org/wpwupingwp/novowrap)
[![PyPI version](https://badge.fury.io/py/novowrap.svg)](https://badge.fury.io/py/novowrap)

# Quick start
Download [the package](https://github.com/wpwupingwp/barcodefinder/releases),
unzip, and run.

__OR__

Open terminal, run
   ```shell
   # Install, using pip (recommended)
   pip install BarcodeFinder --user

   # Initiliaze with Internet
   # Windows
   python -m BarcodeFinder init
   # Linux and MacOS
   python3 -m BarcodeFinder init

   # Run
   # Windows
   python -m BarcodeFinder
   # Linux and MacOS
   python3 -m BarcodeFinder
   ```
# Table of Contents
   * [Quick start](#quickstart)
   * [Feature](#feature)
   * [Prerequisite](#prerequisite)
      * [Hardware](#hardware)
      * [Software](#software)
   * [Installation](#installation)
      * [Portable](#portable)
      * [Install with pip](#Installwithpip)
      * [Install with conda](#Installwithconda)
      * [Initialization](#Initialization)
   * [Usage](#usage)
      * [Quick examples](#quick-examples)
      * [Sequence ID](#sequence-id)
      * [Command line](#commandline)
   * [Input](#input)
   * [Output](#output)
   * [Options](#options)
      * [gb2fasta](#gb2fasta)
      * [evaluate](#evaluate)
      * [primer](#primer)
   * [Performance](#performance)
   * [Citation](#citation)
   * [License](#license)
   * [Q&A](q&a)

# Features
:heavy_check_mark: Automatically collect, organize and clean sequence data from NCBI Genbank
or local: collect data with abundant options; extract CDS,
intergenic spacer, or any other annotations from original sequencep; remove
redundant sequences according to species information; remove invalid or
abnormal sequences/fragments; generate clean dataset with uniform sequence id. 

:heavy_check_mark: Evaluate variance of sequences by calculating nucleotide
diversity, observed resolution, Shannon index, tree resolution, phylogenetic
diversity (original and edited version), gap ratio, and others. Support
sliding-window scanning.

:heavy_check_mark: Design universal primer for the alignment. Support
ambiguous bases in primers.

# Prerequisite
## Hardware
BarcodeFinder requires very few computational resources. A normal PC/laptop is 
enough. For downloading large amount of data, make sure the Internet
connection is stable and fast enough.

## Software
For the portable version, nothing need to be installed manually.

For installing from pip, [Python](https://www.python.org/downloads/) is
required. Notice that the python version should be **3.6** or higher.

:white_check_mark: All third-party dependencies will be automatically
installed with Internet, including `biopython`, `matplotlib`, `coloredlogs`,
`numpy`, `primer3-py`, (python packages), and
[MAFFT](https://mafft.cbrc.jp/alignment/software/),
[IQTREE](http://www.iqtree.org/),
[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download).

# Installation
We assume that users have already installed
[Python3](https://www.python.org/downloads/) (3.7 or above).

## Portable
Download from the [link](https://github.com/wpwupingwp/barcodefinder/releases),
unpack and run with Internet for the first time.
## Install with pip
1. Install [Python](https://www.python.org/downloads/). 3.7 or newer is
   required.
     
2. Open command line, run
```shell
pip install BarcodeFinder --user
```
## Initialization
During the first running, `barcodefinder` will check and initialize the
running environment.  Missing dependencies will be automatically installed.

This step requires Internet connection.
```shell
# Windows
python -m BarcodeFinder init
# Linux and MacOS
python3 -m BarcodeFinder init
```

If BarcodeFinder **FAILED** to install third-party software, please follow these
steps:
For Linux users with root privileges, just use the package manager:
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
For MacOS users with root privileges, install *brew* if it has not been
installed previously:
```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```
If any errors occur, install Xcode from the App Store and retry.

Then:
```
brew install blast mafft brewsci/science/iqtree
```
If using Windows or lacking root privileges, users should follow these
instructions:
1. BLAST+

    * [Windows](https://www.ncbi.nlm.nih.gov/books/NBK52637/)
    * [Linux and MacOS](https://www.ncbi.nlm.nih.gov/books/NBK52640/)
2. MAFFT

    * [Windows](https://mafft.cbrc.jp/alignment/software/windows.html)

        Choose "All-in-one version", download and unzip. Then follow the steps
        in the BLAST+ installation manual to set the _PATH_.
    * [Linux](https://mafft.cbrc.jp/alignment/software/linux.html)

        Choose "Portable package", download and unzip. Then follow the
        instructions of BLAST+ to set the _PATH_ for *MAFFT*.
    * [MacOS](https://mafft.cbrc.jp/alignment/software/macosx.html)

        Choose "All-in-one version", download and unzip. Then follow the steps
        in the BLAST+ installation manual to set the _PATH_.
3. IQ-TREE

    * [Download](http://www.iqtree.org/#download)

        Download the installer according to OS. Unzip and add the path of
        subfolder *bin* to _PATH_.
# Usage
BarcodeFinder is a command-line program. Once a user opens the command line
(Windows) or terminal (Linux and MacOS), just type the command:
```
# Windows
python -m BarcodeFinder [input] -[options] -out [out_folder]
# Linux and MacOS
python3 -m BarcodeFinder [input] -[options] -out [out_folder]
```
## Quick examples
1. Download all _rbcL_ sequences of species in Poaceae family and do
   pre-process.  
```
# Windows
python -m BarcodeFinder.gb2fasta -gene rbcL -taxon Poaceae -out rbcL_Poaceae
# Linux and macOS
python3 -m BarcodeFinder.gb2fasta -gene rbcL -taxon Poaceae -out rbcL_Poaceae
```
2. Download all ITS sequences of _Rosa_ genus. Do pre-process and keep redundant
   sequences:
```
# Windows
python -m BarcodeFinder.gb2fasta -query internal transcribed spacer -taxon Rosa -out Rosa_its -uniq no
# Linux and macOS
python3 -m BarcodeFinder.gb2fasta -query internal transcribed spacer -taxon Rosa -out Rosa_its -uniq no
```
3. Download all Lamiaceae chloroplast genome sequences in the RefSeq database.
   Then do pre-process and evaluation of variance (skip primer designing):
```
# Windows
python -m BarcodeFinder -og cp -refseq -taxon Lamiaceae -out Lamiaceae_cp -skip_primer
# Linux and macOS
python3 -m BarcodeFinder -og cp -refseq -taxon Lamiaceae -out Lamiaceae_cp -skip_primer
```
4. Download sequences of _Zea mays_, set length between 100 bp and 3000 bp,
   and then perform evaluation and primer designing. Note that the space in
   the species name is replaced with underscore "\_".
```
# Windows
python -m BarcodeFinder -taxon Zea_mays -min_len 100 -max_len 3000 -out Zea_mays
# Linux and macOS
python3 -m BarcodeFinder -taxon Zea_mays -min_len 100 -max_len 3000 -out Zea_mays
```
5. Download all _Oryza_ mitochondria genomes, keep the longest sequence for
   each species and run a full analysis: 
```
# Windows
python -m BarcodeFinder -taxon Oryza -og mt -min_len 50000 -max_len 200000 -uniq longest -out Oryza_cp
# Linux and macOS
python3 -m BarcodeFinder -taxon Oryza -og mt -min_len 50000 -max_len 200000 -uniq longest -out Oryza_cp
```
## Input
BarcodeFinder accepts:
1. Genbank queries. Users can use "-query" or combine with other filters;
2. unaligned fasta files. Each file is considered one locus when evaluating
   the variance;
3. alignments (fasta format); and
4. Genbank format files.

## Sequence ID
BarcodeFinder uses a uniform sequence ID for all fasta files that it generates.
```
SeqName|Kingdom|Phylum|Class|Order|Family|Genus|Species|Accession|SpecimenID|Type
# example
rbcL|Poales|Poaceae|Oryza|longistaminata|MF998442|TAN:GB60B-2014|
```
The order of the fields is fixed. The fields are separated by vertical bars
("|"). The space character (" ") was disallowed and was replaced by an
underscore ("\_"). Due to missing data, some fields may be empty. 
* SeqName

    SeqName refers to the name of a sequence. Usually it is the gene name. For
    intergenic spacer, an underscore ("\_") is used to connect two gene's
    names, e.g., "geneA_geneB".

    If a valid sequence name cannot be found in the annotations of the Genbank
    file, BarcodeFinder will use "Unknown" instead.

    For chloroplast genes, if "-rename" option is set, the program will try to
    use regular expressions to fix potential errors in gene names.
* Kingdom

    The kingdom (_Fungi, Viridiplantae, Metazoa_) of a species. For convenience,
    a superkingdom (_Bacteria, Archaea, Eukaryota, Viruses, Viroids_) may be used
    if the kingdom information for a sequence is missing.
* Phylum

    The phylum of the species.
* Class
    
    The class of the species.
    
    Because some species' classes are left emtpy (for instance, basal
    angiosperm) in the cases of plants, BarcodeFinder will guess the class of the
    species.

    Given the taxonomy information in Genbank file:
    ```
    Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
            Spermatophyta; Magnoliophyta; basal Magnoliophyta; Amborellales;
            Amborellaceae; Amborella.
    ```
    BarcodeFinder will use "basal Magnoliophyta" as the class because this
    expression is located before the order ("Amborellales").

* Order

    The order name of the species.
* Family

    The family name of the species.
* Genus

    The genus name of the species, i.e., the first part of the scientific
    name.
* Species

    The specific epithet of the species, i.e., the second part of the
    scientific name of the species. It may contains the subspecies' name.
* Accession

    The Genbank Accession number for the sequence. It does not contain the
    record's version.
* SpecimenID

    The ID of the specimen of the sequence. Usually this value is empty.
* Type

    The type of the sequence. Could be "gene", "spacer", "intron" or else.
## Output
All results will be put in the output folder. If the user does not set the
output path via "-out", BarcodeFinder will create a folder labelled "Result".
* _a_.gb

    The raw Genbank file. The _a_ comes from the query's keyword.
* _a_.plus

    The raw Genbank file plus extended annotations for spacers and introns.
* _a_.fasta

    The converted fasta file of the ".gb" file.
* _b_.primer.csv

    The list of primer pairs in CSV (comma-separated values text) format. The
    _b_ is the name of the locus/fragment (usually a gene or spacer).

    Its title:
    ```
    Locus,Score,Samples,AvgProductLength,StdEV,MinProductLength,MaxProductLength,Coverage,Resolution,TreeValue,AvgTerminalBranchLen,Entropy,LeftSeq,LeftTm,LeftAvgBitscore,LeftAvgMismatch,RightSeq,RightTm,RightAvgBitscore,RightAvgMismatch,DeltaTm,AlnStart,AlnEnd,AvgSeqStart,AvgSeqEnd
    ```

    * Locus

        The name of the locus/fragment.
    * Score

        The score of this pair of primers. Usually the higher, the better.
    * Samples

        The number of sequences which were used to find this pair of primers.

    * AvgProductLength

        The average length of the DNA fragment amplified by this pair of
        primers.
    * StdEV

        The standard deviation of the AvgProductLength. A higher number means
        the primer may amplify different lengths of DNA fragments.
    * MinProductLength

        The minimum length of an amplified fragment.
    * MaxProductLength

        The maximum length of an amplified fragment. Note that all of these
        fields are calculated using given sequences.
    * Coverage

        The coverage of this pair of primers over the sequences it used.
        Calculated with the BLAST result. High coverage means that the pair is
        much more "universal".
    * Resolution

        <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;R_{o}=\frac{n_{uniq}}{n_{total}}" title="R_{o}=\frac{n_{uniq}}{n_{total}}" />

        The *observed resolution* of the sub-alignment sliced by the primer
        pair, which is equal to the number of unique sequences divided by the
        number of total sequences. The value is between 0 and 1.
    * TreeValue

        <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;R_{T}=\frac{n_{internal}}{n_{terminal}}" title="R_{T}=\frac{n_{internal}}{n_{terminal}}" />

        The *tree resolution* of the sub-alignment, which is equal to the
        number of internal nodes on a phylogenetic tree (constructed from the
        alignment) divided by number of terminal nodes. The value is between 0
        and 1.

    * AvgTerminalBranchLen

        The average of the terminal branch's length.
    * Entropy

        <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;E_{H}&space;=&space;\frac{-&space;\sum_{i=1}^{k}{p_{i}&space;\log(p_{i})}}{\log(k)}" title="E_{H} = \frac{- \sum_{i=1}^{k}{p_{i} \log(p_{i})}}{\log(k)}" />

        The Shannon equitability index of the sub-alignment. The value is
        between 0 and 1.
    * LeftSeq

        Sequence of the forward primer. The direction is 5' to 3'.
    * LeftTm

        The melting temperature of the forward primer. The unit is degress
        Celsius (°C).
    * LeftAvgBitscore

        The average raw bitscore of the forward primer, which is calculated by
        BLAST.
    * LeftAvgMismatch

        The average number of mismatched bases of the forward primer, as
        counted by BLAST.
    * RightSeq

        Sequence of reverse primer. The direction is 5' to 3'.
    * RightTm

        The melting temperature of the reverse primer. The unit is degrees
        Celsius (°C).
    * RightAvgBitscore

        The average raw bitscore of the reverse primer, which is calculated by
        BLAST.
    * RightAvgMismatch

        The average number of mismatched bases of the reverse primer, as
        counted by BLAST.
    * DeltaTm

        The difference in the melting temperatures of the forward and reverse
        primers. A pair of primers with a high DeltaTm may result in failure
        during the PCR experiment.
    * AlnStart

        The location of the beginning of the forward primer (5', leftmost of
        primer pairs) in the entire alignment.
    * AlnEnd

        The location of the end of the reverse primer (5', rightmost of primer
        pairs) in the entire alignment.
    * AvgSeqStart

        The average beginning of the forward primer in the original sequences.
        *ONLY USED FOR DEBUG*.
    * AvgSeqEnd

        The average end of the forward primer in the original sequences.
        *ONLY USED FOR DEBUG*.

    The primer pairs are sorted by *Score*. Since the score may not fully
    satisfy the user's specific considerations, it is suggested that primer
    pairs be chosen manually if the first primer pair fails during the PCR
    experiment.
* _b_.primer.fastq

    The fastq format file of a primer's sequence. It contains two sequences,
    and the direction is 5' to 3'. The first is the forward primer, and the
    second is the reverse primer. The quality of each base is equal to its
    proportion of the column in the alignment. Note that the sequence may
    contains amibiguous bases if it was not disabled.
* _b_.pdf

    The PDF format of the figure containing the sliding-window scan result of
    the alignment.
* _b_.png

    The PNG format of the figure containing the sliding-window scan result of
    the alignment.
* _b_.variance.tsv

    The CSV format of the sliding-window scan result. *"Index"* means the
    location of the base in the alignment. Note that the value DOES NOT means
    the variance of the base column; instead it refers to the variance of the
    fragment started from this column.

* Log.txt

    The log file. Contains all the information printed on the screen.
* Options.json

    The JSON file stores all options that the user inputs.
* Loci.csv

    The summary of all loci/fragments, which only contains the variance
    information for each fragment. The only new field, *GapRatio*, means the
    ratio of the gap ("-") in the alignment. A higher value means that the
    sequences may have too many insertions/deletions or the alignment is not
    reliable.
* by_name

    The folder contains *"undivided"* sequences and intermediate results.
    Actually they are "roughly divided" sequences. The original Genbank file
    is firstly divided into different fasta files if the Genbank record
    contains different contents. Usually, one Genbank record contains serveral
    annotated regions (multiple genes, for example). If two records contains
    the same series of annotations (same order), they are put into same fasta
    file. Each file contains the intact sequence form the related Genbank
    record.

* by_gene

    The folder contains divided sequences and intermediate results. After
    the divided step occurred in by_name, BarcodeFinder then divides each
    *cluster* of Genbank records into several fasta files so that each file
    contains only one region (one locus, one gene, one spacer or one
    misc_feature) of the annotation.

    For instance, a record in a "rbcL.gb" file may also contains atpB gene's
    sequences. The "rbcL.fasta" file does not contain any upstream/downstream
    sequences (except for ".expand" files) and "atpB_rbcL.fasta" does not have
    even one base of the atpB or rbcL gene, just the spacer (assuming that the
    annotation is precise).

    User can skip this dividing step by setting "-no_divide" to use the whole
    sequence for analysis. Note that doing so DOES NOT skip the first dividing
    step.

    These two folders can usually be ignored. However, a user may
    utilise one of these intermediate results (especially for those who only
    use BarcodeFinder to collect data from Genbank):
    * _b_.fasta
        The raw fasta file converted directly from the Genbank file containing
        only sequences of one locus/fragment.
    * _b_.expand

        To design primers, BarcodeFinder extend a sequence to its
        upstream/downstream. Users can use "-expand 0" to skip the expansion.
        The next step generates files that all have ".expand" in their
        filenames.
    * _b_.uniq

        Non-redundant sequences.
    * _b_.uniq.aln

        The alignment of the fasta file.
    * _b_.uniq.candidate.fasta

        The candidate primers. This file may contains thousands of records. We do
        not recommend paying attention to it.
    * _b_.uniq.candidate.fastq

        Again, the candidate primers. This time, the file has the quality
        information that equals base's proportion in the column of the
        alignment.
    * _b_.uniq.consensus.fastq

        The fastq format of the consensus sequence of the alignment. Note that
        it contains alignment gap ("-"). Although this may be the most useful
        file in the folder, it is NOT RECOMMENDED that it be used directly
        because consensus-genrating algorithm are optimised for primer design. Hence, the consensus sequence may be different from the _"real"_ consensus.

# Options
## Help
* -h

    Prints help messages for the program. It is highly recommended to use this
    option to see the list of options and their default values.
## General
* -aln filename

    Alignment files that the user provides. The filename can consist of one
    file or a series of files. One can use "?" and "\*" to represent one or
    any characters. *Be sure to use quotation marks*. For example,
    "a\*.alignment" means any file starting with the letter "a" and ending
    with ".alignment".

    It only supports the fasta format. Ambiguous bases and gaps ("-") are
    supported.

* -fasta filename

    User-provided unaligned fasta files. Also supports "\*" and "?". If the
    user wants to use "-uniq" function, the sequences should be renamed. See
    the format for the sequence ID above.

* -gb filename

    User-provided Genbank file or files.

* -stop value

    To stop the running BarcodeFinder at a specific step. BarcodeFinder
    provides an all-in-one solution to find novel DNA barcodes. However, some
    users may only want to use one module. The *value* could be
    * 1
        Only collect data and do pre-processing (download, divide, remove
        redundant, rename, expand); or
    * 2
        Do step 1, and then analyse the variance. Do not design primers.

* -out value

    The output folder's name. All results will be put into the output folder.
    If the user does not set an output path via "-out", BarcodeFinder will
    create a folder named "Result".

    BarcodeFinder does not overwrite the existing folder with the same name.

    It is HIGHLY RECOMMENDED to use only letters, numbers and underscores
    ("_") in the folder name to avoid mysterious errors caused by other
    Unicode characters.

## Genbank
* -allow_mosaic_spacer
    
    If one gene is nested with another gene, normally they do not have spacers.

    However, some users want the fragments between two gene's beginnings and
    ends. This option is for this specific purpose. For normal usage, *do not
    recommend*.
* -allow_repeat

    If genes repeated in downstream, this option will allow the repeat region
    to be extracted, otherwise any repeated region will be omitted.
    
    The default value is False.
* -allow_invert_repeat

    If two genes invert-repeated in downstream, this option will allow the
    spacer of them to be extracted, otherwise the spacer will be omitted.

    For instance, geneA-geneB located in one invert-repeat region (IR) of
    chloroplast genome. In another IR region, there are geneB-geneA. This
    option will extract sequences of two different direction as two unique
    spacers.
    
    The default value is False.
* -email address

    BarcodeFinder uses Biopython to handle the communication between the user
    and the NCBI Genbank database. The database requires that the to provide
    an email address in case of abnormal situations that require NCBI to
    contact the user. The default address was designed to be empty.

    _However, for the convenience of the user, BarcodeFinder will use
    "guest@example.com" if the user does not provide an email address._

* -exclude option

    Use this option to use negative option. For instance, "-exclude Zea
    [organism]" (do not include quotation marks) will add " NOT
    (Zea[organism])" to the query.

    This option can be useful for exclude specific taxon.
    ```
    -taxon Zea -exclude "Zea mays"[organism]
    ```
    This will query all records in genus *Zea* while records of *Zea mays*
    will be exclude.

    For much more complex exclude options, please consider to use "Advance
    search" in Genbank website.
* -gene name

    The gene's name which the user wants to query in Genbank. If the user
    wants to use logical expressions like "OR", "AND", "NOT", s/he should use
    "-query" instead. If there is space in the gene's name, make sure to use
    quotation marks.

    Note that "ITS" is not a gene name--it is "internal transcribed spacer".

    Sometimes "-gene" options may bring in unwanted sequences. For example, if
    a user queries "rbcL[gene]" in Genbank, spacers containing rbcL or rbcL's
    upstream/downstream gene may be found, such as "atpB_rbcL spacer" or atpB.

* -group value

    To restrict a group of species to their *superkingdom* or *kingdom*, the
    value can be

    * animals
    * plants
    * fungi
    * protists
    * bacteria
    * archaea
    * viruses

    It is reported that the "group" filter may return abnormal records, for
    instance, return plants' records when the group is "animal" and the
    "organelle" is "chloroplast". Furthermore, it may match a great number of
    records in Genbank. Hence, we strongly recommend using "-taxon" instead.

    The default *value* is empty.

* -min_len value

    The minimum length of the records downloaded from Genbank. The default
    *value* is 100 (bp). The *number* must be an integer.

* -max_len value

    The maximum length of the records downloaded from Genbank. The default
    *value* is 10000 (bp). The *number* must be an integer.

* -molecular type

    The molecular type, which could be DNA or RNA. The default *type* is empty.

* -og type (or -organelle type)

    Adds "organelle[filter]" to a query to limit results to a given organelle
    type only.

    The *type* could be

    * mitochondrion, or mt
    * plastid, or pl
    * chloroplast, or cp

    Usually, users only want organelle genomes instead of fragments. One
    solution for leaving out fragments is to set "-min_len" and "-max_len" to
    use a length filter to obtain genomes. Another simple solution is to use
    RefSeq only by adding the "-refseq" option.

    For instance,
    ```
    # all chloroplast sequences of Poaceae (not only in RefSeq)
    -taxon Poaceae -og chloroplast -min_len 50000 -max_len 300000
    # all chloroplast sequences of Poaceae (only in RefSeq)
    -taxon Poaceae -og chloroplast -refseq
    ```

    Make sure not to make any typo (e.g., chlorplast, or mitochondria).
* -query string

    The query string provided by the user. It behaves in the same manner as
    the query the user typed into the Search Box in NCBI Genbank's webpage.

    Make sure to follow NCBI's grammar for queries. It can contain several
    words. Remember to add quotation marks if an item contains more than one
    words, for instance, *"Homo sapiens"[organism].

    Do not add quotation marks at the beginning and end of the query string.
    For instance, *"cbs[gene] AND "Homo sapiens"[organism]" may return empty
    results.

* -refseq

    Ask BarcodeFinder to only query sequences in the RefSeq database.
    [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/about/) is considered to be
    of higher quality than the other sequences in Genbank.

    If the user set this option, "-min_seq" and "-max_seq" will be removed to
    set no limit on the sequence length. If the user really wants to set a
    length limit when using "-refseq", s/he can put this filter into the query
    string:
    ```
    # Get all Poaceae sequence in RefSeq, sequences should be 1000 to 10000 bp
    -query refseq[filter] -taxon Poaceae -min_len 1000 -max_len 10000
    ```

    Usually, this option will be combined with "-og" to obtain organelle
    genomes.

    Note that this option is of Boolean type. It IS NOT followed with a *value*.

* -seq_n value

    Download part of records. The *value* should be integer.

    The defaule *value* is None, i.e., download all records.
* -taxon taxonomy

    The taxonomy name. It could be any taxonomic rank from kingdom (same as
    "-group") to species, as long as the user inputs correct name (the
    scientific name of species or taxonomic group in latin, NOT ENGLISH). It
    will restrict the query to the targeted taxonomy unit. Make sure to use
    quotation marks if *taxonomy* has more than one word.

## Pre-process
* -expand value

    The expand length for going upstream/downstream. If set, BarcodeFinder
    will expand the sequence to its upstream/downstream after the dividing
    step to find primer candidates. Set the *number* to 0 to skip.

    The default value is 0 if users set "-stop" to 1 or 2, i.e., users do not
    want to run the primer-design process.

    If users run the whole process but forget to set "-expand", BarcodeFinder
    will automatically set "-expand" to 200. However, users can force the
    program to not to expand the sequence by setting it to 0.
* -max_name_len value

    The maximum length of a feature name. Some annotation's feature name in
    Genbank file is too long, and usually, they are not the target sequence
    the user wants. By setting this option, BarcodeFinder will truncate the
    annotation's feature name if it is too long. By default, the *value* is
    50.
* -max_seq_len value

    The maximum length of a sequence for one annotation. Some annotations'
    sequences are too long (for instance, one gene has two exons, and its
    intron is longer than 10 Kb). This option will skip those long sequences.
    By default, the *value* is 20000 (bp).

    Note that this option is different with "-max_len". This option limits the
    length of one annotation's sequence. The "-max_len" limits the whole
    sequence's length of one Genbank record.

    For an organelle genome's analysis, if the user sets the "-no_divide"
    option, this option will be ignored.
* -no_divide

    If set, it will analyse the whole sequence instead of the divided
    fragments. By default, BarcodeFinder divides one Genbank record into
    several fragments according to its annotation.

    Note that this option is of Boolean type. It IS NOT followed with a *value*.
* -rename

    If set, the program will try to rename genes. For instance, "rbcl" will be
    renamed to "rbcL", and "tRNA UAC" will be renamed to "trnVuac", which
    consists of "trn", the amino acid's letter and transcribed codon. This may
    be helpful if the annotation has nonstandard uppercase/lowercase or naming
    format so it can merge the same sequences to one file for the same locus
    having variant names.

    If using Windows, consider using this option to avoid contradictory
    filenames.

    It is also of Boolean type. The default is not to rename.
* -uniq method

    The method used to remove redundant sequences. BarcodeFinder will remove
    redundant sequences to ensure only one sequence per species by default. A
    user can change its behaviour by setting different methods.
    * longest

        Keep the longest sequence for one species. The program will compare
        the sequence's length from the same species' same locus.
    * random

        The program will randomly pick one sequence for one species' one
        locus if there is more than one sequence.
    * first

        According to the records' order in the original Genbank file, only the
        first sequence of the same species' same locus will be kept. Others
        will be ignored directly. This is the default option due to
        performance considerations.
    * no

        Skip this step. All sequences will be kept.
## Evaluate
* -fast

    If set, BarcodeFinder will skip the calculations for the "tree resolution"
    and "average terminal branch length" to reduce running time.

    Although the program has been optimized greatly, the phylogenetic tree's
    inferences can be time consuming if there are too many species (for
    instance, 10,000). If the user wants to analyse organelle genomes and the
    species are too numerous, the user can set this option to reduce time. The
    "tree resolution" and "average terminal branch length" will both become 0
    in the resultsu file.
* -step value

    The step length for the sliding-window scan. The default *value* is 50. If
    the input dataset is too large, an extreamely small *value* (such as 1 or
    2) may require too much time, especially when the "-fast" option is not
    used.
## Primer design
* -a value

    The maximum number of ambiguous bases allowed in one primer. The default
    *value* is 4.
* -c value

    The minimum coverage of the base and primer. The default *value* is 0.6
    (60%).  It is used to remove primer candidates if its coverage among all
    sequences is smaller than the threshold. The coverage of primers is
    calculated by BLAST.

    Also, it is used to generate a consensus sequence. For one column, if the
    proportion of one type of base (A, T, C, G) is smaller than the threshold,
    the program will try to use an ambiguous base that represents two type of
    bases, and then three, then four ("N").
* -m value

    The maximum number of mismatched bases in a primer. This options is used
    to remove primer candidates if the BLAST results show that there is too
    much mismatch. The default *value* is 4.
* -pmin value

    The minimum length of the primer length. The default *value* is 18.
* -pmax value

    The maximum length of the primer length. The default *value* is 24.
* -r value

    The minimum *observed resolution* of the fragments or primer pairs. The
    default *value* is 0.5. It is used to skip conserved fragments (alignment
    or sub-alignment defined by a pair of primers).

    BarcodeFinder uses the *observed resolution* instead of others for several
    reasons:
    * speed

        The calculation of the *observed resolution* is very fast.
    * accuracy

        Due to the existence of possible alignment errors, the *observed
        resolution* may be higher than the resolutions obtained via other
        evaluation methods. Hence, it is used as a lower bound.  That is to
        say, the program considers that a fragment with a low *observed
        resolution* may not have a satisfactory tree resolution either.

    By setting it to 0, BarcodeFinder can skip this filtration step.
    Meanwhile, the running time may be extremely long.
* -t value

    Only keeps *value* pairs of primers for each highly variant region. The
    default *value* is 1, i.e., only keep the _best_ primer pair. To choose the
    best pairs of primers, the *Score* each pair received is used. To keep
    more pairs, set "-t" to more than 1.
* -tmin value

    The minimum product length (include primer). The default *value* is 300
    (bp). Note this limits the PCR product's length instead of the
    sub-alignment's length.
* -tmax value

    The maximum product length (include primer). The default *value* is 500
    (bp). Note that it limits the length of the PCR product given by the
    primer pair instead of the alignment.

    The "-tmin" and "-tmax" are used to screen primer candidates. It uses
    BLAST results to set the location of primers on each template sequence and
    calculates the average lengths of the products. Because of the variance of
    species, the same locus may have differenct lengths in different species,
    plus with the stretching of the alignment that gaps were added during the
    aligning, please consider adding some *margins* for these two options.

    For instance, if a user wants the amplified length to be smaller than 800
    and greater than 500, s/he could consider setting "-tmin" to 550 and
    "-tmax" to 750.

# Performance
For a taxon that is not very large and includes few fragments, BarcodeFinder
can finish the task in *minutes*. For a large taxon (such as the Asteraceae
family or the whole class of the Poales) and multiple fragments (such as the
chloroplast genomes), the time to complete may be one hour or more on a PC or
laptop.

BarcodeFinder requires less memory (usually less than 0.5 GB, although, for a
large taxon BLAST may require more) and few CPUs (one core is enough). It can
run very well on a normal PC. Multiple CPU cores may be helpful for the
alignment and tree construction steps.

For Windows users, MAFFT [may be very slow due to anti-virus
software](https://mafft.cbrc.jp/alignment/software/windows_without_cygwin.html).
Please consider following [this instruction](https://mafft.cbrc.jp/alignment/software/ubuntu_on_windows.html) to install
Ubuntu on Windows to obtain better results.

# Citation
As yet unpublished.

# License
The software itself is licensed under
[AGPL-3.0](https://github.com/wpwupingwp/novowrap/blob/master/LICENSE) (**not include third-party
software**).

# Q&A
Please submit your questions in the
[Issue](https://github.com/wpwupingwp/barcodefinder/issues) page :smiley:
* Q: I got error message that the program failed to install
  MAFFT/BLAST/IQTREE.

  A: Uncommonly, users in specific area have connection issue for those
  websites. Users have to manually download packages and install (see
  [Software](#software) for the download links).

  For Windows users, please download and unpack files into
  `%HOMEDRIVE%%HOMEPATH%/.barcodefinder`.

  For Linux  and MacOS users, please download and unpack files into
  `~/.barcodefinder`.
