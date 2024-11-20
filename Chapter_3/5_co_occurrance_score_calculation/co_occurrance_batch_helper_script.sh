#!/bin/bash

my_str=$1
my_str=($( basename $my_str ))
my_str=${my_str%_all_global_hits.csv}
global_str=${my_str}_all_global_hits.csv
my_str_local=${my_str}_run_0_139_5_10kb_local_hits.csv
my_str_sequence=${my_str}

# local: sequence_tmp_genome_block_9998.fasta_windows.fastagene_1710499_17_5_10kb_local_hits.csv
# global: sequence_tmp_genome_block_4807.fasta_windows.fastagene_1716578_14_all_global_hits.csv
# sequence: genome_block_9734.fasta_windows.fastagene_1910930_5
python3 co_occurance_calculator.py ${my_str_sequence} local_hits/${my_str_local} global_hits/${global_str} $2