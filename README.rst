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

* Python3 (3.5 or above)
* BLAST+
* IQTREE
* MAFFT

- Biopython
- coloredlogs
- matplotlib
- numpy
- primer3-py

# Project Information

The source code of *BarcodeFinder* is available under AGPLv3 license.
For usage and details, please visit [BarcodeFinder on GitHub](https://github.com/wpwupingwp/BarcodeFinder).
