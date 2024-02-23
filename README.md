[![PyPI version](https://badge.fury.io/py/OGU.svg)](https://badge.fury.io/py/OGU)

# Quick start

- Install Python 3 (3.9 or newer).
- Open terminal, run

   ```shell
   # Install, using pip (recommended)
   pip3 install OGU --user

   # Initialize with Internet
   # Windows
   python -m OGU init
   # Linux and macOS
   python3 -m OGU init

   # Run
   # Windows
   python -m OGU
   # Linux and macOS
   python3 -m OGU
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

:heavy_check_mark: Automatically collect, organize and clean sequence data
from NCBI GenBank or local: collect data with abundant options; extract CDS,
intergenic spacer, or any other annotations from original sequence; remove
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

`Organelle Genome Utilities (OGU)` requires very few computational resources.
A normal PC/laptop is enough. For downloading large amount of data, make sure
the Internet connection is stable and fast enough.

## Software

For the portable version, nothing need to be installed manually.

For installing from pip, [Python](https://www.python.org/downloads/) is
required. Notice that the python version should be higher than **3.6**.

:white_check_mark: All third-party dependencies will be automatically
installed with Internet, including `biopython`, `matplotlib`, `coloredlogs`,
`numpy`, `primer3-py`, (python packages), and
[MAFFT](https://mafft.cbrc.jp/alignment/software/),
[IQTREE](http://www.iqtree.org/),
[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download).

# Installation

We assume that users have already installed
[Python3](https://www.python.org/downloads/) (3.9 or above).

## Install with pip

1. Install [Python](https://www.python.org/downloads/). 3.9 or newer is
   required.

2. Open command line, run

```shell
pip3 install OGU --user
```

## Initialization

During the first running, `OGU` will check and initialize the
running environment. Missing dependencies will be automatically installed.

This step requires Internet connection.

By default, the program will automatically finish initialization, if any error
occurs, users can choose one of the following methods:

### Automatic

Run the following command.

```shell
# Windows
python -m OGU init
# Linux and macOS
python3 -m OGU init
```

### Use prepared package
According to your system, download related compressed file from [packages](https://github.com/wpwupingwp/OGU/releases).

For Windows users: 
1. paste `%HOMEDRIVE%%HOMEPATH%/` to explorer's address bar and open.
2. create `.OGU` folder. Don't miss the dot.
3. open `.OGU` folder, paste downloaded compressed file and unzip. Make sure after
decompress there are three folders in `.OGU`.

For Linux and macOS users, please download and unpack files into
`~/.OGU`.

### Manually install

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

For macOS users with root privileges, install `brew` if it has not been
installed previously:

```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

If any errors occur, install `Xcode-select` and retry.

Then:

```
brew install blast mafft brewsci/science/iqtree
```

If using Windows or lacking root privileges, users should follow these
instructions:

1. BLAST+

    * [Windows](https://www.ncbi.nlm.nih.gov/books/NBK52637/)
    * [Linux and macOS](https://www.ncbi.nlm.nih.gov/books/NBK52640/)
2. MAFFT

    * [Windows](https://mafft.cbrc.jp/alignment/software/windows.html)

      Choose "All-in-one version", download and unzip. Then follow the steps
      in the BLAST+ installation manual to set the `PATH`.
    * [Linux](https://mafft.cbrc.jp/alignment/software/linux.html)

      Choose "Portable package", download and unzip. Then follow the
      instructions of BLAST+ to set the `PATH` for `MAFFT`.
    * [macOS](https://mafft.cbrc.jp/alignment/software/macosx.html)

      Choose "All-in-one version", download and unzip. Then follow the steps
      in the BLAST+ installation manual to set the `PATH`.
3. IQ-TREE

    * [Download](http://www.iqtree.org/#download)

      Download the installer according to OS. Unzip and add the path of
      subfolder `bin` to `PATH`.

# Usage

## Graphical user interface
Open the command line (Windows) or terminal (Linux and macOS),
run

```bash
OGU
```

or 
```bash
# linux and macos
python3 -m OGU
# windows
python -m OGU
```

## command line
Once a user opens the command line (Windows) or terminal (Linux and macOS), 
just type the command:

```
# Windows
python -m OGU [input] -[options] -out [out_folder]
# Linux and macOS
python3 -m OGU [input] -[options] -out [out_folder]
```

## Quick examples

1. Download all `rbcL` sequences of species in Poaceae family and do
   pre-process.

```
# Windows
python -m OGU.gb2fasta -gene rbcL -taxon Poaceae -out rbcL_Poaceae
# Linux and macOS
python3 -m OGU.gb2fasta -gene rbcL -taxon Poaceae -out rbcL_Poaceae
```

2. Download all ITS sequences of _Rosa_ genus. Do pre-process and keep redundant
   sequences:

```
# Windows
python -m OGU.gb2fasta -query internal transcribed spacer -taxon Rosa -out Rosa_its -uniq no
# Linux and macOS
python3 -m OGU.gb2fasta -query internal transcribed spacer -taxon Rosa -out Rosa_its -uniq no
```

3. Download all Lamiaceae chloroplast genomic sequences in the RefSeq database.
   Then do pre-process and evaluation of variance (skip primer designing):

```
# Windows
python -m OGU -og cp -refseq yes -taxon Lamiaceae -out Lamiaceae_cp
# Linux and macOS
python3 -m OGU -og cp -refseq yes -taxon Lamiaceae -out Lamiaceae_cp
```

4. Download sequences of _Zea mays_, set length between 100 bp and 3000 bp,
   and then perform evaluation and primer designing. Note that the space in
   the species name is replaced with underscore "\_".

```
# Windows
python -m OGU -taxon Zea_mays -min_len 100 -max_len 3000 -out Zea_mays -primer
# Linux and macOS
python3 -m OGU -taxon Zea_mays -min_len 100 -max_len 3000 -out Zea_mays -primer
```

5. Download all _Oryza_ mitochondria genomes in RefSeq database, keep the
   longest sequence for each species and run a full analysis:

```
# Windows
python -m OGU -taxon Oryza -og mt -min_len 50000 -max_len 200000 -uniq longest -out Oryza_cp -refseq yes -primer
# Linux and macOS
python3 -m OGU -taxon Oryza -og mt -min_len 50000 -max_len 200000 -uniq longest -out Oryza_cp -refseq yes -primer
```

## Sequence ID

`Organelle Genome Utilities` uses a uniform sequence id format for input fasta files and all output
sequences.

```
Locus|Kingdom|Phylum|Class|Order|Family|Genus|Species|Accession|SpecimenID_Isolate|Type
# example
rbcL|Viridiplantae|Streptophyta|Magnoliopsida|Poales|Poaceae|Oryza|longistaminata|MF998442|TAN:GB60B-2014|gene
```

The order of the fields is fixed. The fields are separated by vertical bars
("|"). The space character (" ") was disallowed and was replaced by an
underscore ("\_"). Due to missing data, some fields may be empty.

`Locus`: SeqName refers to the name of a sequence. Usually it is the gene
name. For intergenic spacer, an underscore ("\_") is used to connect two
gene's names, e.g., "geneA_geneB".

If a valid sequence name cannot be found in the annotations of the GenBank
file, `Organelle Genome Utilities` will use "Unknown" instead.

For chloroplast genes, if "-rename" option is set, the program will try to use
regular expressions to fix potential errors in gene names.

`Kingdom`: The kingdom (_Fungi, Viridiplantae, Metazoa_) of a species. For
convenience, a superkingdom (_Bacteria, Archaea, Eukaryota, Viruses, Viroids_)
may be used if the kingdom information for a sequence is missing.

`Phylum`: The phylum of the species.

`Class`: The class of the species.

Because some species' classes are empty (for instance, basal angiosperm),
`Organelle Genome Utilities` will guess the class of the species.

Given the taxonomy information in GenBank file:

```
Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
    Spermatophyta; Magnoliophyta; basal Magnoliophyta; Amborellales;
    Amborellaceae; Amborella.
```

`Organelle Genome Utilities` will use "basal Magnoliophyta" as the class because this
expression locates before the order name ("Amborellales").

`Order`: The order name of the species.

`Family`: The family name of the species.

`Genus`: The genus name of the species, i.e., the first part of the scientific
name.

`Species`: The specific epithet of the species, i.e., the second part of the
scientific name of the species. It may contain the subspecies' name.

`Accession`: The GenBank Accession number for the sequence. It does not
contain the record's version.

`SpecimenID` and `Isolate`: Specimen ID and Isolate ID of the sequence. May be empty.

`Type`: Type of the sequence. It is usually "gene" or "spacer".

## Command line

:exclamation: In Linux and macOS, Python2 is `python2` and Python3 is
`python3`. However, in Windows, Python3 is called `python`, too. Please
notice the difference.

* Show help information of each module

 ```shell
 # Windows
 python -m OGU -h
 python -m OGU.gb2fasta -h
 python -m OGU.evaluate -h
 python -m OGU.primer -h
 # Linux and macOS
 python3 -m OGU.gb2fasta -h
 python3 -m OGU.evaluate -h
 python3 -m OGU.primer -h
 ```

* Full process

 ```shell
 # Windows
 python -m OGU -gene [gene name] -taxon [taxon name] -og [organelle type] -out [output name]
 # Linux and macOS
 python3 -m OGU -gene [gene name] -taxon [taxon name] -og [organelle type] -out [output name]
 ```

* Collect, convert, and clean GenBank data with gb2fasta module

 ```shell
 # Windows
 python -m OGU.GB2fasta -gene [gene name] -taxon [taxon name] -og [organelle type] -out [output name]
 # Linux and macOS
 python3 -m OGU.gb2fasta -gene [gene name] -taxon [taxon name] -og [organelle type] -out [output name]
 ```

* Evaluate variance of given fasta files

 ```shell
 # Windows
 python -m OGU.evaluate -fasta [fasta files]
 # Linux and macOS
 python3 -m OGU.evaluate -fasta [input file]
 ```

* Design universal primers of given alignments.

 ```shell
 # Windows
 python -m OGU.primer -aln [alignment files]
 # Linux and macOS
 python3 -m OGU.primer -aln [alignment files]
 ```

## Visualize
 
Here are to jupyter notebooks for visualize analyze result as detailed circle
figure:
  - `Visualize/draw_circle_plastid.ipynb`: for plastid genomes
  - `Visualize/draw_circle_mitochondria.ipynb`: for mitochondria genomes

Since users may want to customize the figure, we provide jupyter notebooks 
instead of packaged code. Users can get the result following these steps.
0. Run `pip install jupyterlab` to install jupyter notebooks
1. Double click to open in jupyter notebook, Visual Studio Code or other IDEs you prefer.
2. Edit `filename` to the Evaluation.csv you got
3. Edit `gb_file` to extended gb file you got. Remember to generate it with 
`-out_debug` in OGU.gb2fasta
4. If visualize plastid data, you need provide `LSC, SSC, IRa, IRb` lengths. Or
you can use default value, which is for *Tobacum*.
5. Edit color themes as your wish
6. Run all cells to output pdf figure
 
# Input

`Organelle Genome Utilities` accepts:

1. GenBank queries. Users can use "-query" or combine with any other filters;
2. GenBank-format files.
3. Unaligned fasta files. Each file is considered as one locus when evaluating
   the variance;
4. Alignments (fasta format).

# Output

All results will be put in the output folder. If the user does not set the
output path via "-out", `Organelle Genome Utilities` will create a folder labelled "Result".

In the output folder, several sub-folders will be created.

* GenBank

  Raw GenBank files.

* Divide

  Fasta files converted from the GenBank file. Each file represents a
  fragment of the original sequence according to the annotation.

  For instance, a record in a "rbcL.gb" file may also contain atpB gene's
  sequences. The "rbcL.fasta" file does not contain any upstream/downstream
  sequences and "atpB_rbcL.fasta" does not have even one base of the atpB or
  rbcL gene, just the spacer (assuming the annotation is precise).

  User can skip this dividing step with the option "-no_divide".
* Fasta

  Raw fasta files users provided.
* Unique

  Fasta files after removing redundant sequences.
* Expanded_fasta

  To design primers, `Organelle Genome Utilities` extend a sequence to its
  upstream/downstream. Only used in the primer module.
* Alignment

  Aligned fasta files.

  `.aln`: The aligned fasta files.

  `.-consensus.fastq`: The fastq format of the consensus sequence of the
  alignment. Note that it contains alignment gap ("-"). It is NOT
  RECOMMENDED to be used directly because the consensus-generating algorithm is
  optimised for primer design.
* Evaluate

  Including output files from the evaluation module.

  `.pdf`: The PDF format of the figure containing the sliding-window scan
  result of the alignment.

  `.csv`: The CSV format file of the sliding-window scan result. `"Index"`
  means the location of the base in the alignment.

* Primer

  Including output files from the primer module.

  `.primer.fastq`: The fastq format file of a primer's sequence. It contains
  two sequences, and the direction is 5' to 3'. The first is the forward
  primer, and the second is the reverse primer. The quality of each base is
  equal to its proportion of the column in the alignment. Note that the
  sequence may contain ambiguous bases if it was not disabled.

  `.primers.csv`: The list of primer pairs in CSV (comma-separated values
  text) format.

  `.candidate.fasta`: The candidate primers. This file may contain
  thousands of records. Do not recommend paying attention to it.

  `.candidate.fastq`: Again, the candidate primers. This time, each file has
  the quality information that equals to the proportion of the bases in the
  column of the alignment.

* Temp

  Including temporary files. Could be safely deleted .

In the output folder, there are some other important output files:

* Primers.csv

  The list of primer pairs in CSV (comma-separated values text) format.

  Its title:
    ```
    Locus,Samples,Score,AvgProductLength,StdEV,MinProductLength,MaxProductLength,Coverage,Observed_Res,Tree_Res,PD_terminal,Entropy,LeftSeq,LeftTm,LeftAvgBitscore,LeftAvgMismatch,RightSeq,RightTm,RightAvgBitscore,RightAvgMismatch,DeltaTm,AlnStart,AlnEnd,AvgSeqStart,AvgSeqEnd
    ```

  `Locus`: The name of the locus/fragment.

  `Samples`: The number of sequences used to find this pair of primers.

  `Score`: The score of this pair of primers. Usually the higher, the better.

  `AvgProductLength`: The average length of the DNA fragment amplified by
  this pair of primers.

  `StdEV`: The standard deviation of the AvgProductLength. A higher number
  means the primer may amplify different lengths of DNA fragments.

  `MinProductLength`: The minimum length of an amplified fragment.

  `MaxProductLength`: The maximum length of an amplified fragment. Note that
  all of these fields are calculated using given sequences.

  `Coverage`: The coverage of this pair of primers over the sequences it
  used. Calculated with the BLAST result. High coverage means that the pair
  is much more "universal".

  `Observed_Res`: The `observed resolution` of the sub-alignment sliced by
  the primer pair, which is equal to the number of unique sequences divided
  by the number of total sequences. The value is between 0 and 1.

    <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;R_{o}=\frac{n_{uniq}}{n_{total}}" title="R_{o}=\frac{n_{uniq}}{n_{total}}" />

  `Tree_Res`: The `tree resolution` of the sub-alignment, which is equal to
  the number of internal nodes on a phylogenetic tree (constructed from the
  alignment) divided by number of terminal nodes. The value is between 0 and
    1.

    <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;R_{T}=\frac{n_{internal}}{n_{terminal}}" title="R_{T}=\frac{n_{internal}}{n_{terminal}}" />

  `PD_terminal`: The average of the terminal branch's length. It's an edited
  version of the `Phylogenetic Diversity` for DNA barcoding evaluation.

  `Entropy`: The `Shannon equitability` index of the sub-alignment. The value
  is between 0 and 1.

    <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;E_{H}&space;=&space;\frac{-&space;\sum_{i=1}^{k}{p_{i}&space;\log(p_{i})}}{\log(k)}" title="E_{H} = \frac{- \sum_{i=1}^{k}{p_{i} \log(p_{i})}}{\log(k)}" />

  `LeftSeq`: Sequence of the forward primer. The direction is 5' to 3'.

  `LeftTm`: The melting temperature of the forward primer. The unit is
  degree Celsius (°C).

  `LeftAvgBitscore`: The average raw bitscore of the forward primer, which
  is calculated by BLAST.

  `LeftAvgMismatch`: The average number of mismatched bases of the forward
  primer, as counted by BLAST.

  `RightSeq`: Sequence of reverse primer. The direction is 5' to 3'.

  `RightTm`: The melting temperature of the reverse primer. The unit is
  degrees Celsius (°C).

  `RightAvgBitscore`: The average raw bitscore of the reverse primer, which
  is calculated by BLAST.

  `RightAvgMismatch`: The average number of mismatched bases of the reverse
  primer, as counted by BLAST.

  `DeltaTm`: The difference in the melting temperatures of the forward and
  reverse primers. A pair of primers with a high DeltaTm may result in
  failure during the PCR experiment.

  `AlnStart`: The location of the beginning of the forward primer (5',
  leftmost of primer pairs) in the entire alignment.

  `AlnEnd`: The location of the end of the reverse primer (5', rightmost of
  primer pairs) in the entire alignment.

  `AvgSeqStart`: The average beginning of the forward primer in the original
  sequences.  *ONLY USED FOR DEBUG*.

  `AvgSeqEnd`: The average end of the forward primer in the original
  sequences.  *ONLY USED FOR DEBUG*.

  The primer pairs are sorted by `Score`. Since the score may not fully
  satisfy the user's specific considerations, it is suggested that primer
  pairs be chosen manually if the first primer pair fails during the PCR
  experiment.

* Log.txt

  The log file. Contains all the information printed on the screen.

* Evaluation.csv

  The summary of all loci/fragments, which only contains the variance
  information for each fragment. One of the new field, `GapRatio`, means the
  ratio of the gap ("-") in the alignment. `PD` means the original version
  of the phylogenetic diversity and `PD_stem` means an alternative version
  of it which only calculate the length of the stem branch in the
  phylogenetic tree.

# Options

Here are some general options for the program and submodule:

`-h`: Prints help messages of the program or one of the module.

`-gb [filename]`: User-provided GenBank file or files. Could be one or more
files that separated by space.

For instance,

```
# one file
-gb sequence.gb
# multiple files
-gb matK.gb rbcL.gb Oryza.gb Homo_sapiens.gb
```

`-fasta [filename]`: User-provided unaligned fasta files. Could be one or
multiple.

`-aln [filename]`: Alignment files that the user provides. Could be one or
multiple.

It only supports the fasta format. Ambiguous bases and gaps ("-") are supported.

`-out [folder name]`: The output folder's name. All results will be put into
the output folder. If the user does not set an output path via "-out",
`Organelle Genome Utilities` will create a folder named "Result".

`OGU` does not overwrite the existing folder with the same name.

It is HIGHLY RECOMMENDED to use only letters, numbers and underscores ("\_") in
the folder name to avoid mysterious errors caused by other Unicode characters.

Options below are for specific modules.

## gb2fasta

### Query

Options used for querying NCBI GenBank.

`-taxon [taxonomy name]`: The taxonomy name. It could be any taxonomic rank
from kingdom (same as "-group") to species, as long as the user inputs correct
name (the scientific name of species or taxonomic group in latin, NOT
ENGLISH). It will restrict the query to the targeted taxonomy unit. Make sure
to use quotation marks if `taxonomy` has more than one word or use underscore
to replace space, for instance `"Zea mays"` or `Zea_mays`.

`-gene [gene name]`: The gene's name which the user wants to query in GenBank.
If the user wants to use logical expressions like "OR", "AND", "NOT", s/he
should use "-query" instead. If there is space in the gene's name, make sure
to use quotation marks.

Note that "ITS" is not a gene name--it is "internal transcribed spacer".

Sometimes "-gene" options may bring in unwanted sequences. For example, if a
user queries "rbcL[gene]" in GenBank, spacer sequences may contain _rbcL_ or
_rbcL_'s upstream/downstream gene, such as "atpB_rbcL spacer" or _atpB_.

`-og [ignore|both|no|mt|mitochondrion|cp|chloroplast|pl|plastid]`: Query
organelle sequences or not. The default value is `ignore`.

    - `ignore`: do not consider organelle type, same as GenBank website's
      default setting.

    - `both`: only query organelle sequences, including both plastid and
      mitochondrion.

    - `no`: exclude organelle sequences from the query.

    - `cp` or `chloroplast` or `pl` or `plastid`: only query plastid sequences

    - `mt` or `mitochondrion`: only query mitochondrion sequences.

`-refseq [both|yes|no]`: query in RefSeq database or not. The default value is
`both`.

    - `both`: query all sequences in or not in RefSeq database, same as NCBI
      website's default setting.

    - `yes`: only query sequences in RefSeq database.

    - `no`: exclude sequences in RefSeq database.

[RefSeq](https://www.ncbi.nlm.nih.gov/refseq/about/) is considered to have
higher sequence and annotation quality than GenBank. This option could be used
for getting nuclear/organelle genomes from NCBI. In this situation (`-refseq
yes`), the length limit will be removed automatically.

`-count [number]`: Restrict numbers of sequences to be downloaded. The default
value `0` means no restriction.

`-min_len [length]`: The minimum length of the records downloaded from
GenBank. The default value is `100` (bp). The number must be an integer.

`-max_len [length]`: The maximum length of the records downloaded from
GenBank. The default value is `10000` (bp). The number must be an integer.

`-date_start [yyyy/mm/dd]`: The beginning of the release data range of the
sequences, the format is yyyy/mm/dd.

`-date_end [yyyy/mm/dd]`: The end of the release data range of the sequences,
the format is yyyy/mm/dd.

`-molecular [all|DNA|RNA]`: The molecular type,
which could be DNA or RNA. The
default is `all`--no restriction.

`-email [email address]`: NCBI GenBank database requires users to provide
an email address in case of abnormal situations that NCBI need to contact
the user. For convenience, `OGU` will use
"guest@example.com" if the user does not provide an email address. _However_,
it is better to provide a real email address for potential contact.

`-query [expression]`: The query string provided by the user. It behaves in
the same manner as the query the user typed into the Search Box in NCBI
GenBank's webpage.

Make sure to follow NCBI's grammar for queries. It can contain several words.
Remember to add quotation marks if an item contains more than one words, for
instance, `"Homo sapiens"[organism]`, or use underscore to replace space,
`Homo_sapiens[organism]`.

`-exclude [expression]`: Use this option to use negative option. For instance,
"-exclude Zea [organism]" (do not include quotation marks) will add " NOT
(Zea[organism])" to the query.

This option can be useful for excluding a specific taxon.

```
-taxon Zea -exclude "Zea mays"[organism]
```

This will query all records in genus *Zea* while records of *Zea mays* will be
excluded.

For much more complex exclude options, please consider to use "Advance search"
in GenBank website.

`-group [all|animals|plants|fungi|protists|bacteria|archaea|viruses]`: To
restrict the query in given group. The default value is `all`--no
restriction.

It is reported that the "group" filter may return abnormal records, for
instance, return plants' records when the group is "animal" and the
"organelle" is "chloroplast". Furthermore, it may match a great number of
records in GenBank. Hence, we strongly recommend using "-taxon" instead.

### Divide

Options used for converting GenBank files to fasta files.

`-out_debug`: If you are going to use visualize pipeline to draw detailed circle
figure, use this option to generate extended version genbank file.

`-no_divide`: If set, it will analyse the whole sequence instead of the
divided fragments. By default, `OGU` divides one GenBank record into
several fragments according to its annotation.

`-rename`: If set, the program will try to rename genes. For instance, "rbcl"
will be renamed to "rbcL", and "tRNA UAC" will be renamed to "trnVuac", which
consists of "trn", the amino acid's letter and transcribed codon. This may be
helpful if the annotation has nonstandard uppercase/lowercase or naming format
that it can merge the same sequences to one file for the same locus having
variant names.

If using Windows operating system, consider using this option to avoid
contradictory filenames.

`-unique [longest|first|no]`: The method used to remove redundant sequences.
`OGU` will remove redundant sequences to ensure only one sequence per
species by default. A user can change its behaviour by setting different
methods.

    - `first`: According to the records' order in the original GenBank file,
      only the first sequence of the same species' same locus will be kept.
      Others will be ignored directly. This is the default option due to
      performance considerations.

    - `longest`: Keep the longest sequence for one species. The program will
      compare the sequence's length from the same species' same locus.

    - `no`: Skip this step. All sequences will be kept.

`-allow_mosaic_spacer`: If one gene nested with another gene, normally they
do not have spacers. The default value is `False`.

However, some users want the fragments between two gene's beginnings and ends.
This option is for this specific purpose (e.g., matK-trnK_UUU). For normal
usage, *do not recommend*.

`-expand [number]`: The expansion length in upstream/downstream. If set,
`OGU` will expand the sequence to its upstream/downstream after the
dividing step to find primer candidates. The default value is `0`.

Note that this option is different with "-max_len". This option limits the
length of one annotation's sequence. The "-max_len" limits the whole
sequence's length of one GenBank record.

`-allow_repeat`: If genes repeated in downstream, this option will allow the
repeat region to be extracted, otherwise any repeated region will be omitted.
The default value is `False`.

`-allow_invert_repeat`: If two genes invert-repeated in downstream, this
option will allow spacers of them to be extracted, otherwise the spacer
will be omitted. The default value is `False`.

For instance, geneA-geneB located in one invert-repeat region (IR) of
chloroplast genome. In another IR region, there are geneB-geneA. This option
will extract sequences of two different direction as two unique spacers.

`-max_name_len [number]`: The maximum length of a feature name. Some
annotation's feature name in GenBank file is too long, and usually, they are
not the target sequence the user wants. By setting this option, `OGU`
will truncate the annotation's feature name if it is too long. By default, the
value is `50`.

`-max_gene_len [value]`: The maximum length of a sequence for one annotation.
Some annotations' sequences are too long (for instance, one gene has two
exons, and its intron is longer than 10 Kb). This option will skip those long
sequences. By default, the value is `20000` (bp).

## Evaluate

`-ig` or `-ignore_gap`: ignore gaps in the alignment.

`-iab` or `-ignore_ambigous`: ignore ambiguous bases in the alignment.

`-quick`: skip sliding-window scan.

`-size [number]`: the window size of the sliding window scan. The default
value is `500`.

`-step [number]`: the step size of the sliding window scan. The default value
is `50`.

`-skip_primer`: skip primer designing. The default value is `False`.

## Primer design

`-coverage [value]`: The minimum coverage of the base and primer. The default
value is `0.5` (50%). It is used to remove primer candidates if its coverage
among all sequences is smaller than the threshold. The coverage of primers is
calculated by BLAST.

`-res [value]`: The minimum *observed resolution* of the fragments or primer
pairs. The default *value* is 0.3 (30%). The value should be in 0.0 to 1.0.

`OGU` uses the *observed resolution* instead of others because of the
speed. Also, it is considered to be the lower bound of the real resolution
that a fragment with a low *observed resolution* may not have a satisfactory
tree resolution/phylogenetic diversity, either.

`-pmin [length]`: The minimal length of the primer. The default *value* is 20.

`-pmax [length]`: The maximal length of the primer. The default *value* is 25.

`-topn [number]`: How many pairs of primers is kept for each input alignment.
The default value is `1`, i.e., only keep the _best_ primer pair according to
its `score`. To keep more pairs, set "-t" to more than 1.

`-amin [length]`: The minimum amplified length (include primer). The default
value is `300` (bp). Note this limits the PCR product's length instead of the
sub-alignment's length.

`-amax [length]`: The maximum amplified length (include primer). The default
value is `800` (bp).

The "-amin" and "-amax" are used to screen primer candidates. It uses BLAST
results to set the location of primers on each template sequence and
calculates the average lengths of the products. Because of the variance of
species, the same locus may have different lengths in different species, plus
with the stretching of the alignment that gaps were added during the aligning,
please consider adding some *margins* for these two options.

For instance, if a user wants the amplified length to be smaller than 800 and
greater than 500, s/he could consider setting "-amin" to 550 and "-amax" to

750.

`-ambiguous [number]`: The maximum number of ambiguous bases allowed in one
primer. The default value is `4`.

`-mismatch [number]`: The maximum number of mismatched bases in a primer. This
option is used to remove primer candidates if the BLAST results show that
there is too much mismatch. The default value is `4`.

# Performance

For a taxon that is not very large and includes few fragments, `OGU`
can finish the task in *minutes*. For a large taxon (such as the Asteraceae
family or the whole class of the Poales) containing multiple fragments (such
as the chloroplast genomes), the time to complete may be one hour or more on a
PC or laptop.

`OGU` requires few memories (usually less than 0.5 GB, although, for a
large taxon BLAST may require more) and few CPUs (one core is enough). It can
run very well on a normal PC. Multiple CPU cores may be helpful for the
alignment and tree construction steps.

For Windows users, MAFFT [may be very slow due to antivirus
software](https://mafft.cbrc.jp/alignment/software/windows_without_cygwin.html).
Please consider
following [this instruction](https://mafft.cbrc.jp/alignment/software/ubuntu_on_windows.html) to
install
Ubuntu on Windows to obtain better results.

# Citation

As yet unpublished.

# License

The software itself is licensed under
[AGPL-3.0](https://github.com/wpwupingwp/OGU/blob/master/LICENSE) (**not include
third-party
software**).

# Q&A

Please submit your questions in the
[Issue](https://github.com/wpwupingwp/OGU/issues) page :smiley:
* Q: The first-time run is slow.

  A: OGU will automaticlly install third-party software (MAFFT/BLAST/IQTREE)
  from AWS at first-time running. Microsoft Windows users, especially in some 
  regions may have slow connection. Please be patient, or you can manually 
  install them. See [Initialization](#Initialization).