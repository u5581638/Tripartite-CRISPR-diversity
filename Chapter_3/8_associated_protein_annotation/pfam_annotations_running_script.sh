#!/bin/bash

#PBS -P do77
#PBS -q normalsr
#PBS -l walltime=8:00:00
#PBS -l mem=120GB
#PBS -l ncpus=48
#PBS -l wd
#PBS -l storage=scratch/xi88+gdata/xi88+gdata/va71

module load python3/3.10.4

# script to use HMMscan/HHblits to search for remote homology to each representative protein sequence in the DEFLOC and Pfam databases respectively

# Pfam run
find split_sequences/ -name "*" -type "f" | xargs -n 1 -I {} -P 48 /g/data/va71/alphafold/install/hh-suite/build/src/hhblits -i {} -o {}_pfam_descriptions.txt  -cpu 1 -e 0.001 -maxseq 10 -d /g/data/va71/alphafold/datasets/alphafold2/pfam2/pfam -v 0

find split_sequences/ -name "*" -type "f" | xargs -n 1 -I {} -P 48 python3 custom_hhsuite_parser.py {}_pfam_descriptions.txt

# DEFLOC run
find split_sequences/ -name "*" -type "f" | xargs -n 1 -I {} -P 48 /g/data/va71/alphafold/install/hmmer-3.3.2/src/hmmscan -o {}_crisprcas_hits.txt --cpu 1 /g/data/va71/alphafold/datasets/alphafold2/crisprcasfinder_hmm_profiles/non_dup_merged_all_models_master_TIGR.HMM {}

find split_sequences/ -name "*" -type "f" | xargs -n 1 -I {} -P 48 python3 custom_hmm_scan_parser.py {}_crisprcas_hits.txt

