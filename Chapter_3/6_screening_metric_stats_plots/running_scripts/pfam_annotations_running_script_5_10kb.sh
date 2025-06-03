#!/bin/bash

#PBS -P do77
#PBS -q normalsr
#PBS -l walltime=8:00:00
#PBS -l mem=120GB
#PBS -l ncpus=48
#PBS -l wd
#PBS -l storage=scratch/xi88+gdata/xi88+gdata/va71

module load python3/3.10.4
python3 sequence_splitter_5_10kb.py all_run_5_10kb_sequences.fasta

find split_sequences_5_10kb/ -name "*" -type "f" | xargs -n 1 -I {} -P 48 /g/data/va71/alphafold/install/hh-suite/build/src/hhblits -i {} -o {}_pfam_descriptions.txt  -cpu 1 -e 0.001 -maxseq 10 -d /g/data/va71/alphafold/datasets/alphafold2/pfam2/pfam -v 0

find split_sequences_5_10kb/ -name "*" -type "f" | xargs -n 1 -I {} -P 48 python3 custom_hhsuite_parser.py {}_pfam_descriptions.txt

# hhsuite_parser.pdb70_parse_pfam(folder_name + "hmmer_sequence_files/" + sequence_id + "_hhblits_hits.txt.fa")

find split_sequences_5_10kb/ -name "*" -type "f" | xargs -n 1 -I {} -P 48 /g/data/va71/alphafold/install/hmmer-3.3.2/src/hmmscan -o {}_crisprcas_hits.txt --cpu 1 /g/data/va71/alphafold/datasets/alphafold2/crisprcasfinder_hmm_profiles/non_dup_merged_all_models_master_TIGR.HMM {}

find split_sequences_5_10kb/ -name "*" -type "f" | xargs -n 1 -I {} -P 48 python3 custom_hmm_scan_parser.py {}_crisprcas_hits.txt

# phmmer_frames = hmm_scan_parser.hmmscan_parse(folder_name + "hmmer_sequence_files/" + sequence_id + "_hmmscan_hits.fa")