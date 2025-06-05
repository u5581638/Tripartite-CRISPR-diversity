#!/bin/bash
#PBS -P va71
#PBS -q normalsr
#PBS -l walltime=48:00:00
#PBS -l mem=200GB
#PBS -l ncpus=100
#PBS -l wd
#PBS -l storage=scratch/xi88+gdata/xi88+gdata/va71

# script used to run IQtree2
# INPUT: Multiple sequence alignment of representative cas10d protein orthologs (after prior clustering and extracting representatives using mmseqs2)
# OUTPUT: An approximate maximum-likelihood tree generated via IQtree2
export PATH=$PATH:/g/data/va71/gaetan_data/iqtree-2.2.2.6-Linux/bin

my_input=typeID_cas10d.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta_rep_aa_clustered.fa_aln.fasta

iqtree2 -s $my_input -st AA -m MFP -T 100 -alrt 1000 -B 1000 -bnni

