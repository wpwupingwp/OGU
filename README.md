# About
Here are some python scripts I wrote. Most of them process fasta/fastq/gb
files.

From now on, I will add some descriptions for each program.

You can type 
>python3 pyfile -h

or

>python3 pyfile

to print usage of each program.

You can ask me any question about these programs via
**wpwupingwp@outlook.com** .
# Requirement
* python3
* Bioython

And notice that all scripts were just tested on Linux system, although it may works fine on Windows.

## filter.py

Use info in **wanted** to filter fasta file, the output sequences contain **wanted** in their id.

## get_whole.py

Give fragment sequence, return whole length of specific input fragment from given big fasta file via BLAST.

## rename.py

Generate fasta file with special id format from given genbank format file ID looks like this:

>organism|specimen_voucher|isolate|gene_name_series

## vlookup_assistant.py

Expand a given table according to range.

Input table (CSV format) looks like this:

>    A,B,C

It will generate a new table:

>    D,E 

where D was expanded from range(B, C) and E is related A.
