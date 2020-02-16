#!/usr/bin/env python3


#                   USAGE:
# Puth the script file in the same directory with the data table file
# Don't forget to put arguments after script name:
# first argument - data table file name (e.g. 'Supplementary data')
# second argument - prefix in the output files (e.g. '16-feb-2020_nt_seq_list')
# third argument (optional) is the output directory,
# which is the folder name that should exist in the current directory (e.g. 'nt_seq_list_by_gene')


import pandas as pd
import numpy as np
import re
import sys
import os



if len(sys.argv) <= 2:
    print("ERROR: don't forget to type input file name as the first argument")
    print("ERROR: don't forget to type desired output prefix as the second argument")
    exit(1)
    quit()

input_dataframe = sys.argv[1]
prefix = sys.argv[2]
out_dir = ""
out_dir = sys.argv[3]

filename_re = '.+' + prefix + ".fasta"
for filename in os.listdir("./"):
    if re.search(filename_re, filename) is not None:
        print("Error: file with this prefix already exists, try another prefix")
        quit()

df = pd.read_csv(input_dataframe, sep='\t', encoding='UTF-8')
#df = pd.read_csv(input_dataframe)
#print(df)

colnames = df.columns.values.tolist()

genome_parts = list()

for colname in colnames:
    gnm_prt = re.findall('(.+)_start', colname)
    if len(gnm_prt) != 0:
        genome_parts.append(gnm_prt[0])
    continue
#print(genome_parts)


utrs = ["5-UTR","3-UTR"]


for index, row in df.iterrows():
    if len(row["Sequence_identifier"]) == 0:
        break

    fasta_short_header = ">" + row["Sequence_identifier"]+ "\n"
    fasta_short_header = fasta_short_header.replace(" ", "_")

    # get full genome sequence first
    full_genome = row["Sequence"]
    full_genome = full_genome.upper()
    full_genome = full_genome.replace("U", "T")

    #create fasta file
    genome_fasta = fasta_short_header + full_genome + "\n"
    fasta_genome_name =  "./" + out_dir + "/" + "genome" + "_" + prefix + ".fasta"
    fasta_genome = open(fasta_genome_name, "a+")
    fasta_genome.write(genome_fasta)
    # get ORF
    orf_start = row["C_start"]
    orf_stop = row["NS5_stop"]
    if not np.isnan(orf_start).any() and not np.isnan(orf_stop).any():
        orf_start_coord = int(orf_start) - 1
        orf_stop_coord = int(orf_stop)
        orf = full_genome[orf_start_coord:orf_stop_coord]
        #verify sequence length
        if len(orf)%3 != 0:
            print("!!! TO CUT FROM 3'-UTR --- ", len(orf)%3, " nt")
            orf = orf[:(len(orf)-len(orf)%3)]
        orf_fasta = fasta_short_header + orf + "\n"
        print("Virus", row["Sequence_identifier"], "***\n")
        print("Full genome is ", len(full_genome), " nt long, full ORF is ",
         len(orf), " nt long", "\n")

        #create fasta file
        fasta_orf_name =  "./" + out_dir + "/" + "orf" + "_" + prefix + ".fasta"
        fasta_orf = open(fasta_orf_name, "a+")
        fasta_orf.write(orf_fasta)

    #get separate genes
    for g in genome_parts:
        start_colname = g + "_start"
        start = row[start_colname]
        stop_colname = g + "_stop"
        stop = row[stop_colname]
        if not np.isnan(start).any() and not np.isnan(stop).any():
            start = int(start)
            stop = int(stop)
            start_coord = start - 1
            stop_coord = stop

            # make sure the gene starts with at least from 1
            if start_coord >= 0 and stop_coord >= 0:
                full_gene = full_genome[start_coord:stop_coord]
                #verify sequence length
                if g not in utrs and len(full_gene)%3 != 0:
                    print("!!! TO CUT FROM 3'-UTR --- ", len(full_gene)%3, " nt")
                    full_gene = full_gene[:(len(full_gene)-len(full_gene)%3)]
                gene_fasta = fasta_short_header + full_gene + "\n"
                fasta_gene_name = "./" + out_dir + "/" + g + "_" + prefix + ".fasta"
                fasta_gene = open(fasta_gene_name, "a+")
                fasta_gene.write(gene_fasta)
