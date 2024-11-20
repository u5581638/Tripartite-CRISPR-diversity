#!/bin/bash

#PBS -P va71
#PBS -q normalsr
#PBS -l walltime=24:00:00
#PBS -l mem=350GB
#PBS -l ncpus=4
#PBS -l wd
#PBS -l storage=scratch/xi88+gdata/xi88+gdata/va71

module load python3/3.10.4

find global_hits/ -name "sequence_tmp_genome_block*" -type "f" | xargs -n 1 -I {} -P 1 ./co_occurrance_batch_helper_script.sh {} all_5_10kb_run_0_139_co_occurrance_scores_table.csv