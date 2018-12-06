# Table of Contents
   * [Introduction](#introduction)
   * [Prerequisite](#prerequisite)
      * [Software](#software)
      * [Python module](#python-module)
      * [Internet](#internet)
   * [Installation](#installation)
   * [Usage](#usage)
      * [Quick examples](#quick-examples)
      * [Input](#input)
      * [Sequence ID](#sequence-id)
      * [Output](#output)
   * [Options](#options)
      * [Help](#help)
      * [General](#general)
      * [Genbank](#genbank)
      * [Preprocess](#preprocess)
      * [Evaluate](#evaluate)
      * [Primer Design](#primer-design)
   * [Performance](#performance)
# Introduction
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
    Pi, Shannon Index, observed resolution, tree resolution and average
    terminal branch length, etc. If the result is lower than given threshold,
    i.e., it does not have efficient resolution, this alignment were skipped.

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
# Prerequisite
## Software
* Python3 (3.5 or above)
* BLAST+
* IQTREE
* MAFFT
## Python module
* Biopython
* matplotlib
* numpy
* primer3-py
## Internet
The data retrieve function requires Internet connection. Please ensure you
have stable network and reasonable Internet traffic charge for downloading
large size of data.
# Installation
Assume that you already installed [Python3](https://www.python.org/downloads/)
(3.5 or above, *3.6* is recommended). Firstly, install BarcodeFinder.

The easiest way is to use pip. Make sure your pip is not out of date (18.0 or
newer), then
```
# As administator
pip3 install BarcodeFinder
# Normal user
pip3 install BarcodeFinder --user
```
For some version of python (eg. python 3.7 or above), pip may ask user to
provide compiler. If you do not have (especially for Windows user), we
recommend [this](https://www.lfd.uci.edu/~gohlke/pythonlibs/) website. You can
download compiled wheel file and use pip to install it:
```
# As administator
pip3 install wheel_file_name
# Normal user
pip3 install wheel_file_name --user
```
Secondly, you need to install dependent software. BarcodeFinder has assistant
function to automatically install dependent software if it cannot find them,
i.e., you can skip this step as you wish. However, it is *highly recommended*
to follow official installation procedure to make it easy for management and
give you a clean working directory.

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
For MacOS user, assume you have root privileges, if you do not have *brew*,
install it:
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

        Choose "All-in-one version", download and unzip. Then follow the step
        in BLAST+ installation manual to set _PATH_.
    * [Linux](https://mafft.cbrc.jp/alignment/software/linux.html)

        Choose "Portable package", download and unzip. Then follow the
        instruction of BLAST+ to set _PATH_.
    * [MacOS](https://mafft.cbrc.jp/alignment/software/macosx.html)

        Choose "All-in-one version", download and unzip. Then follow the step
        in BLAST+ installation manual to set _PATH_.
3. IQTREE

    * [Download](http://www.iqtree.org/#download)

        Download installer according to your OS. Unzip and add the path of
        subfolder *bin* into _PATH_
# Usage
BarcodeFinder is a command line program. Once you open your command line
(Windows) or terminal (Linux and MacOS), just type the command:
```
# Windows
python -m BarcodeFinder [input] -[options] -out [out_folder]
# Linux and MacOS
python3 -m BarcodeFinder [input] -[options] -out [out_folder]
```
## Quick examples
1. Download all _rbcL_ sequences of plants and do pre-process:
```
# Windows
python -m BarcodeFinder -query rbcL -group plants -stop 1 -out rbcL
# Linux and macOS
python3 -m BarcodeFinder -query rbcL -group plants -stop 1 -out rbcL
```
2. Download all ITS sequences of _Rosa_ and do pre-process:
```
# Windows
python -m BarcodeFinder -query "internal transcribed spacer" -taxon Rosa -stop 1 -out Rosa_its
# Linux and macOS
python3 -m BarcodeFinder -query "internal transcribed spacer" -taxon Rosa -stop 1 -out Rosa_its
```
3. Download all Rosaceae chloroplast genome sequences in RefSeq database, plus
   your own data.  Then do pre-process and evaluation of variance (do not
   design primers):
```
# Windows
python -m BarcodeFinder -organelle chloroplast -refseq -taxon Rosaceae -out Poaceae_cpg -fasta my_data.fasta -stop 2
# Linux and macOS
python3 -m BarcodeFinder -organelle chloroplast -refseq -taxon Rosaceae -out Poaceae_cpg -fasta my_data.fasta stop 2
```
4. Download sequences of _Zea mays_, set length between 100 bp and 3000 bp,
   plus your aligned data, then do full analysis:
```
# Windows
python -m BarcodeFinder -taxon "Zea mays" -min_len 100 -max_len 3000 -out Zea_mays -aln my_data.aln
# Linux and macOS
python3 -m BarcodeFinder -taxon "Zea mays" -min_len 100 -max_len 3000 -out Zea_mays -aln my_data.aln
```
5. Download all _Oryza_ chloroplast genomes (not only in RefSeq database),
   keep the longest sequence for each species and do full analysis:
```
# Windows
python -m BarcodeFinder -taxon Oryza -organelle chloroplast -min_len 50000 -max_len 500000 -uniq longest -out Oryza_cp
# Linux and macOS
python3 -m BarcodeFinder -taxon Oryza -organelle chloroplast -min_len 50000 -max_len 500000 -uniq longest -out Oryza_cp
```
## Input
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
## Sequence ID
BarcodeFinder use uniform sequence ID in all fasta files it generated.
```
SeqName|Order|Family|Genus|Species|Accession|SpecimenID
# example
rbcL|Poales|Poaceae|Oryza|longistaminata|MF998442|TAN:GB60B-2014
```
The order of the seven fields is fixed. Each field was seperated by the
vertical bar ("|"). The space character (" ") was disallowed and it was
replaced by underscore ("\_").  Because of technical issue, the order rank may
be empty for animals.
* SeqName

    It means the name of sequences. Usually it its the gene name. For spacer,
    it is "geneA_geneB" that use underscore ("\_") to connect two gene's name.

    Note that the original gene name in genbank file were renamed to try to
    fix part of annotation error.

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
## Output
All results will be put in the output folder. If you didn't set output path by
"-out", BarcodeFinder will create a folder named by current time, for example,
"2018-11-19T16-41-59.330217".
* _a_.gb

    The raw genbank file. The _a_ comes from the keyword of query.
* _a_.fasta

    The converted fasta file of the ".gb" file.
* _b_.primer.csv

    The list of primer pairs. CSV (comma-separated values
    text) format. The _b_ is the name of the locus/fragment (usually gene or
    spacer).

    Its title:
    ```
    Locus,Score,Sequences,AvgProductLength,StdEV,MinProductLength,MaxProductLength,Coverage,Resolution,TreeValue,AvgTerminalBranchLen,Entropy,LeftSeq,LeftTm,LeftAvgBitscore,LeftAvgMismatch,RightSeq,RightTm,RightAvgBitscore,RightAvgMismatch,DeltaTm,AlnStart,AlnEnd,AvgSeqStart,AvgSeqEnd
    ```

    * Locus

        The name of locus/fragment.
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
* _b_.primer.fastq

    The fastq format file of primer's sequence. It contains two sequences and
    the direction is 5' to 3'. The first is the forward primer and the second
    is the reverse primer. The quality of each base equal to its proportion of
    the column in the alignment. Note that it may contains amibiguous base if
    you did not disable it.
* _b_.pdf

    The PDF format of the figure of sliding-window scan result of the
    alignment.
* _b_.png

    The PNG format of the figure of sliding-window scan result of the
    alignment.
* _b_.variance.tsv

    The CSV format of sliding-window scan result. The *"Index"* means the
    location of the base in the alignment. Note that the value DO NOT means
    the variance of the column of the base but the fragment started from this
    column.

* Log.txt

    The log file. Contains all information printed in screen.
* Options.json

    The JSON file stored all options user inputed.
* Loci.csv

    The summary of all loci/fragments which only contains the variance
    information of each fragment. The only new field, *GapRatio*, means the
    ratio of the gap ("-") in the alignment. Higher value means the sequences
    may be too variaty that the alignment is not reliable.
* by_name

    The folder contains *"undivided"* sequences and intermediate results.
    Actually it is "roughly divided" sequences. The original genbank file was
    firstly divided into different fasta files if the genbank record contained
    different content. Usualy one genbank record contains serveral annotated
    region (multiple genes, for example). If two records contains same series
    of annotation (same order), they were put into same fasta file. Each file
    contains the intact sequence in the related genbank record.

* by_gene

    The folder contains divided sequences and intermediate results. After
    divided step occured in by_name, BarcodeFinder then divided each
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
        The raw fasta file directly converted from genbank file which only
        contains sequences of one locus/fragment.
    * _b_.expand

        To design primers, BarcodeFinder automaticly extend the sequences to
        its upstream/downstream. Users can use "-expand 0" to skip the
        expansion. The following step generates files all have ".expand" in
        their filename.
    * _b_.uniq

        Redundant sequences were removed.
    * _b_.uniq.aln

        The alignment of the fasta file.
    * _b_.uniq.candidate.fasta

        The candidate primers. It may contains thousands of records. Do not
        recommend to pay attention to it.
    * _b_.uniq.candidate.fastq

        Still, the candidate primers. This time it has quality which equals to
        base's proportion in the column of the alignment.
    * _b_.uniq.consensus.fastq

        The fastq format of the consensus sequence of the alignment.  Note
        that it contains aligment gap ("-").  Although this may be the most
        useful file in the folder, it is NOT RECOMMENDED to directly use it as
        the consensus sequence of the alignment because the
        consensus-genrating algorithm were optimized for primer design that it
        may be different with the _"real"_ consensus.

# Options
## Help
* -h

    Print help message of the program. It is highly recommended to use this
    option to see the usage of options and their default value.
## General
* -aln filename

    Alignment files user provided. The filename could be one file, or a series
    of files. You can use "?" and "\*" to represent one or any characters. *Be
    sure to use quotation mark* to quote it. For example, "a\*.alignment"
    means any file start with letter "a" and end with ".alignment".

    Only support fasta format. Ambiguous base and gap ("-") were supported.

* -fasta filename

    User provided unaligned fasta files. Also support "\*" and "?". If you
    want to use "-uniq" function, please rename your sequences.  See the
    format of sequence ID above.

* -gb filename

    User provided genbank file or files.

* -stop value

    To break the running of BarcodeFinder in the specific step. BarcodeFinder
    provided all-in-one solution to find novel DNA barcode.  However, some
    user may only want to use one module. The *value* could be
    * 1
        Only collect data and do preprocess (download, divide, remove
        redundant, rename, expand).
    * 2
        Do step 1, and then analyze the variance. Do not design primers.
    * 3
        Do all things (which have no effect).

* -out value

    The output folder's name. All results will be put in the output folder. If
    you didn't set output path by "-out", BarcodeFinder will create a folder
    named by current time, for example, "2018-11-19T16-41-59.330217".

    BarcodeFinder does not overwrite existing folder with same name.

    It is HIGHLY RECOMMENDED to use only letter, number and underscore ("_")
    in the folder name to avoid mysterious error caused by other Unicode
    characters. :)

## Genbank
* -email address

    BarcodeFinder use Biopython to handle the communication between user and
    NCBI Genbank database. The database requires user to provide an email
    address in case of abnormal situation that NCBI want to contact you. The
    default address is empty.

    _However, for convenience of the user, BarcodeFinder will use
    "guest@example.com" if user did not provide the email._

* -gene name

    The gene's name user wants to query in Genbank. The "OR" and "AND" were
    allowed. Make sure to use quotation mark. For instance, "atpB OR rbcL"
    (include quotation mark) means gene atpB or rbcL.

    HOWEVER, the "-gene" options usually does not behave like what we thought.
    If you query "rbcL[gene]" in Genbank, most of the time, NCBI still gives
    you other gene's records.

* -group value

    To restrict group of species, or the *kingdom*, the value could be

    * animals
    * plants
    * fungi
    * protists
    * bacteria
    * archaea
    * viruses

    Make sure do not have typo.

    If use this option, BarcodeFinder may possibly download a huge number of
    records from Genbank. Please ensure you have stable and unexpensive
    Internet connection. The default *value* is empty.

* -min_len value

    The minium length of the records downloaded from Genbank.  The default
    *value* is 100 (bp). The *number* must be integer.

* -max_len value

    The maximum length of the records downloaded from Genbank.  The default
    *value* is 10000 (bp). The *number* must be integer.

* -molecular type

    The molecular type, could be DNA or RNA. The default *type* is empty.

* -organelle type

    Add "organelle[filter]" in query that sequence was limited on given
    organelle type. The *type* could be

    * mitochondrion
    * plastid
    * chloroplast

    Usually users only want organelle genome instead of fragment. One solution
    is to set "-min_len" and "-max_len" to use length filter to get genomes.
    Another simple solution is to only use RefSeq by adding "-refseq" option.

    For instance,
    ```
    # all chloroplast sequences of Poaceae (not only in RefSeq)
    -taxon Poaceae -organelle chloroplast -min_len 50000 -max_len 300000
    # all chloroplast sequences of Poaceae (only in RefSeq)
    -taxon Poaceae -organelle chloroplast -refseq
    ```

    Make sure do not have typo.
* -query string

    The query string user provied. It behaves same with the query you typed in
    the Search Box in NCBI Genbank's webpage.

    Make sure to follow NCBI's grammer of query. Please do not forget
    quotation mark if it has more than one word. The default *string* is
    empty.

* -refseq

    Ask BarcodeFinder to only query sequences in RefSeq database.
    [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/about/) was
    considered to have higher quality than other sequences in Genbank.

    If use this option, "-min_seq" and "-max_seq" will be removed to set no
    limit on sequence length. If user really want to set length limit when
    using "-refseq", user could put this filter into query string:
    ```
    # Get all Poaceae sequence in RefSeq, sequences should be 1000 to 10000 bp
    -query refseq[filter] -taxon Poaceae -min_len 1000 -max_len 10000
    ```

    Usually, this option will be combined with "-organelle" to get high quality organelle genomes.

    Note that this option is boolen type. It DOES NOT followed with a *value*.
* -taxon taxonomy

    The taxonomy name. It could be any taxonomy rank. From kingdom (same with
    "-group") to species or subspecies, as long as you input correct name
    (scientific name of species or taxonomic group, latin, NOT ENGLISH), it
    will restriced query to your target taxonomy unit. Make sure to use
    quotation mark if *taxonomy* has more than one word.

## Preprocess
* -expand value

    The expand length of upstream/downstream. The default *value* is 200 (bp).
    If set, BarcodeFinder will expand the sequence to its upstream/downstream
    after dividing step to find primer candidates. Set the *number* to 0 to
    skip.
* -max_name_len value

    The maximum length of feature name. Some annotation's feature name  in
    genbank file is too long and usually they are not target sequence user
    wanted. By setting this option, Barcodefinder will truncate annotation's
    feature name if too long. By default the *value* is 50.
* -max_seq_len value

    The maximum length of sequence of one annotation. Some annotation's
    sequence is too long (for instance, one gene have two exons and its intron
    is longer than 10 Kb), This option will skip those long sequences. By
    default the *value* is 20000 (20 KB).

    Note that this option is different with "-max_len". This option limits the
    length of one annotation's sequence. The "-max_len" limits the whole
    sequence's length of one genbank record.

    For organelle genome's analysis, if you set "-no_divide" option, this
    option will be ignored.
* -no_divide

    If set, analyze whole sequence instead of divided fragment. By default,
    BarcodeFinder divided one genbank records to several fragments according
    to its annotation.

    Note that this option is boolen type. It DOES NOT followed with a *value*.
* -rename

    If set, the program will try to rename genes. For instance, "rbcl" will be
    renamed to "rbcL", and "tRNA UAC" will be renamed to "trnVuac", which
    consists of "trn", the amino acid's letter and transcribed codon. This may
    be helpful if the annotation has nonstandard uppercase/lowercase or naming
    format. So it can merge same sequences to one file which are same locus
    but have variant name.

    If you use Windows, please use this option to avoid confliction of
    filename.

    It is also a boolen type. The default is not to rename.
* -uniq method

    The method to remove redundant sequences. BarcodeFinder will remove
    redundant sequences to ensure only one sequence for one species by
    default. You can change its behaviour by set different method.
    * longest

        Keep the longest sequence for one species. The program will compare
        sequence's length for the same species' same locus.
    * random

        The program will randomly pick one sequence for one species's one
        locus if there are more than one sequence.
    * first

        According to records' order in the original genbank file, only the
        first sequence of the same species' same locus will be kept. Others
        will be directly ignored. This is the default option because of the
        consideration of performance.
    * no

        Skip this step, all sequences will be kept.
## Evaluate
* -f

    If set, BarcodeFinder will skip the calculation of "tree resolution" and
    "average terminal branch length" to reduce running time.

    Although the program has been optimized greatly, the phylogenetic tree's
    inference could be time-consuming if there are too many species (for
    instance, 10000). If you want to analyze organelle genome and the species
    are too many, you can set this option to reduce time. The "tree
    resolution" and "average terminal branch length" will become 0 in the
    result file.
* -s value

    The step length for sliding-window scan. The default *value* is 50. If
    the input data is too big, extreamly small *value* (such as 1 or 2) may
    cause too much time, especially when the "-f" option were not used at the
    same time.
## Primer Design
* -a value

    The maximum number of ambiguous bases allowed in one primer. The default
    *value* is 4.
* -c value

    The minimum coverage of base and primer. The default *value* is 0.6 (60%).
    It was used to remove primer candidates if its coverage among all
    sequences were smaller than threshold. The coverage of primers were
    calculated by BLAST.

    Also it was used to generate consensus sequence. For one column, if the
    proportion of one type of base (A, T, C, G) is smaller than the threshold,
    the program will try to use ambiguous base which represent two type of
    bases, and then three, then four ("N").
* -m value

    The maximum number of mismatch bases in primer. This options was used to
    remove primer candidates if the BLAST results show that it has too much
    mismatch. The default *value* is 4.
* -pmin value

    The minimum length of the primer length. The default *value* is 18.
* -pmax value

    The maximum length of the primer length. The default *value* is 24.
* -r value

    The minimum *observed resolution* of the fragments or primer pairs. The
    default *value* is 0.5. It was used to skip conserved fragment (alignment
    or sub-alignment defined by a pair of primer).

    BarcodeFinder use *observed resolution* instead of others is because,
    * speed

        The calculation of *observed resolution* is very fast.
    * accuracy

        Because of the exist of possible alignment error, the *observed
        resolution* may be higher than other evaluate method's result. Hence
        it was used to be the lower bound. That is to say, the program
        considers that if a fragment has low *observed resolution*, its *tree
        resolution* may not satisfy the requirement, either.

    By set it to 0, BarcodeFinder can skip this step of filteration.
    Meanwhile, the running time may be high.
* -t value

    Only keep *value* pairs of primers for each high varient region. The
    default *value* is 1, i.e., only keep the _best_ primer pair. Which is the
    best is according to *Score* it got as desciribed before. To get much more
    choice you can set "-t" to more than 1.
* -tmin value

    The minimum product length (include primer). The default *value* is 300
    (bp). Note this limit the PCR product's length instead of sub-alignment's
    length.
* -tmax value

    The maximum product length (include primer). The default *value* is 500
    (bp). Note that it limit the length of PCR product given by the primer
    pair instead of the alignment.

    The "-tmin" and "-tmax" were used to screen primer candidates. It use
    BLAST result to set the location of primer on each template sequence and
    calculate average length of the product. Because of the variance
    that same locus may have differenct length in different species, plus
    with the strecthing of the alignment that gaps were added during the
    aligning, please consider to add some *margin* for these two options.

    For instance, if you want the amplified length smaller than 800 and
    greater than 500, maybe you could consider to set "-tmin" to 550 and
    "-tmax" to 750.

# Performance
For taxon that not very large and few fragments, BarcodeFinder could finish
the task in *minutes*. For large taxon (such as Asteraceae family or the whole
plants kingdom) and multiple fragments (such as chloroplast genomes) the time
may be one hour or more -- in a PC/laptop.

BarcodeFinder requires few memory (usually less than 1GB, for large taxon
BLAST may require more) and CPU (one core is enough). It can run in normal PC
very well. Multiple CPU cores may be helpful for the alignment and tree
construction steps.

For Windows user, it's said that MAFFT [may be very slow due to anti-virus
software](https://mafft.cbrc.jp/alignment/software/windows_without_cygwin.html).
Please consider to follow [this instruction] (https://mafft.cbrc.jp/alignment/software/ubuntu_on_windows.html) to install Ubuntu on Windows and get better experiment.
