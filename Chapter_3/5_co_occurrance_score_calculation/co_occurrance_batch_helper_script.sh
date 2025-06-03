#!/bin/bash

# helper script to compute co_occurrance scores
my_str=$1
my_str=($( basename $my_str ))
my_str=${my_str%_all_global_hits.csv}
global_str=${my_str}_all_global_hits.csv
my_str_local=${my_str}_run_0_139_5_10kb_local_hits.csv
my_str_sequence=${my_str}

# INPUT: local and global decoded protein sequence hits from tblastn (in csv format)
# i.e. 1. local: sequence_tmp_genome_block_9998.fasta_windows.fastagene_1710499_17_5_10kb_local_hits.csv
# 	   2. global: sequence_tmp_genome_block_4807.fasta_windows.fastagene_1716578_14_all_global_hits.csv
# OUTPUT: N/A (parameters passed to co_occurrance_calculator.py)
# SHELL: run through co_occurrance_calculation_script.sh

python3 co_occurance_calculator.py ${my_str_sequence} local_hits/${my_str_local} global_hits/${global_str} $2