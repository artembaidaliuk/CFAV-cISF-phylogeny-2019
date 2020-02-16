# CFAV-cISF-phylogeny-2019

This repository contains nucleotide sequence allignments and metadata for sequences used to reconstruct the phylogeny of cell-fusing fusing agent virus (CFAV) and classical insect-specific flaviviruses (cISFs).

The alignment files correspond to three data subsets for CFAV phylogeny and one for cISF:

* Data subset i - whole genome sequences of CFAV. Alignment files contain 'fgo' in the names.
* Data subset ii - reduced data subset with selected whole genome sequences of CFAV for tree topology hypothesis testing. Alignment files contain 'red' in the names.
* Data subset iii - whole genome and partial genome sequences of CFAV. Alignment files contain 'pge' in the names.
* Data subset iv - whole genome sequences of cISFs. 

Partition files are added where gene partition model were tested.

A python script is also included in order to split the sequences into separate gene (or UTR) sequence lists directly from the data table. The sequence lists could then be used to perform the alignments from scratch. In order to create sequence lists for different data subsets make sure to subset the data table before running the script.

Example:

$ python3 genesplit_custom_v0.0.1.py Supplementary_data.txt 16-feb-2020 nt_seq_list_by_gene

where arg 1 is the data table file in the same directory as the script, arg 2 is an arbitrary prefix for the output sequence list files (in this case just the date), arg 3 is the name of the output folder (should already exist in the current directory).

