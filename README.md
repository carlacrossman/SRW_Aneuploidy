## Detection of sex chromosomal aneuploidy (XXY) in a Southern Right Whale (Eubalaena australis) using read depth and coverage from whole genome sequence data

This repository contains data and scripts associated with *Detection of sex chromosomal aneuploidy (XXY) in a Southern Right Whale (Eubalaena australis) using read depth and coverage from whole genome sequence data*. It contains the following:

- **ZFX-ZFY_fasta** : This directory contains a fasta file of the ZFX-ZFY amplicons for each sample used in this study 
- **Aneuploidy_plots.Rmd** :  This R markdown notebook contains the code used to generate the figures.
- **merge_mark_coverage.sh** : This is the second script to run that merges bam files, marks duplicate reads and calculates coverage and depth across each chromosome.
- **prep_fastq_2_bam.sh** : This is the first script to run that trims and maps the fastq data to the reference genome

The raw data used in this study are archived in NCBI’s Sequence Read Archive under BioProject PRJNA914998.