#!/bin/bash
#PBS -P va71
#PBS -q normalsr
#PBS -l walltime=48:00:00
#PBS -l mem=200GB
#PBS -l ncpus=104
#PBS -l wd
#PBS -l storage=scratch/xi88+gdata/xi88+gdata/va71

# script used to run IQtree2

export PATH=$PATH:/g/data/va71/gaetan_data/iqtree-2.2.2.6-Linux/bin

my_input=cas12a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta_rep_aa_clustered.fa_aln.fasta

iqtree2 -s $my_input -st AA -m MFP -T 104 -alrt 1000 -B 1000 -bnni

