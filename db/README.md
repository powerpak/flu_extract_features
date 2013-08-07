# README: db folder

This folder is for files related to the database--schemas for setting it up, data files to be loaded into it, etc.

Some of these files should be retrieved from IRD, and can be updated as needed.  They are as follows:

## `ird_strains.tsv`

Created by going to [IRD Strain Search][ird_ss], clicking Search with default settings, selecting all ~200k segments, clicking Download, and downloading as a tab delimited file.  It will download with an .xls extension, but it is actually a TSV file.

[ird_ss]: http://www.fludb.org/brc/influenza_sequence_search_strain_display.do?method=ShowCleanSearch&decorator=influenza

## Use Rake to generate these

1. `PROTEIN-A/`
2. `protein-a.fasta`
3. `blastdb/`