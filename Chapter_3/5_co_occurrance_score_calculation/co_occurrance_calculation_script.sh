#!/bin/bash

#PBS -P va71
#PBS -q normalsr
#PBS -l walltime=24:00:00
#PBS -l mem=350GB
#PBS -l ncpus=4
#PBS -l wd
#PBS -l storage=scratch/xi88+gdata/xi88+gdata/va71

module load python3/3.10.4

# compute the co_occurrance scores based on the number of local (hits within 20kb of a CRISPR-array) hits / global hits in all 9.8TB of sequencing data
# INPUT: individual protein sequence hits (both local and global). These must use a related file naming convention for both to be parsed to co_occurrance_batch_helper_script.sh and onward to co-occurrance_calculator.py
# i.e. sequence_tmp_genome_block_1.fasta_global_hits.csv
# OUTPUT: table of co-occurrance scores
# i.e. all_5_10kb_run_0_139_co_occurrance_scores_table.csv (or as specified at input)
# SHELL: 
find global_hits/ -name "sequence_tmp_genome_block*" -type "f" | xargs -n 1 -I {} -P 1 ./co_occurrance_batch_helper_script.sh {} all_5_10kb_run_0_139_co_occurrance_scores_table.csv