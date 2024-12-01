#!/bin/bash

#PBS -P va71
#PBS -q normalsr
#PBS -l walltime=12:00:00
#PBS -l mem=250GB
#PBS -l ncpus=8
#PBS -l wd
#PBS -l storage=scratch/xi88+gdata/xi88+gdata/va71


module load python3/3.10.4
cd run_0_10_local_window_search/
# decode tBLASTn results from a search against DNA within 20kb of predicted CRISPR-arrays.
python3 ../sequence_decoder.py 0_10_5_10kb_local_hits.csv
