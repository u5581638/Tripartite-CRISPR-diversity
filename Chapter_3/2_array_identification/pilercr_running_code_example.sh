#!/bin/bash
#PBS -P do77
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -l mem=96GB
#PBS -l ncpus=48
#PBS -l wd
#PBS -l storage=scratch/xi88+gdata/xi88+gdata/va71

# Code to identify CRISPR-arrays in genome blocks using PILER-CR
# Note: This is a reproduction of the original code which was run interactively via command line
# INPUT: labelled genome blocks (by block name) in FASTA format.
# i.e. genome_block_1.fasta_labelled.fasta (the suffix may have been removed)
# OUTPUT: PILERCR array position file in quasi-FASTA format (as.txt)
# SHELL: (see below)

find labelled_genomes/ -name "genome_block_*" | xargs -n 1 -I {IN} -P 1 /g/data/va71/crispr_pipeline_annotation/pilercr1.06/pilercr -in {IN} -out {IN}.txt -quiet