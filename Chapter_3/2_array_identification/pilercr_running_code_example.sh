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

find labelled_genomes/ -name "genome_block_*" | xargs -n 1 -I {} -P 1 /g/data/va71/crispr_pipeline_annotation/pilercr1.06/pilercr -in {} -out {}.txt -quiet