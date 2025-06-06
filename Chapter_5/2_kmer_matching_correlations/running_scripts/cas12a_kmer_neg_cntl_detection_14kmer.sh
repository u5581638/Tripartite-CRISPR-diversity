#!/bin/bash
#PBS -P do77
#PBS -q normalsr
#PBS -l walltime=10:00:00
#PBS -l mem=250GB
#PBS -l ncpus=52
#PBS -l wd
#PBS -l storage=scratch/xi88+gdata/xi88+gdata/va71

module load python3/3.10.4
module load blast
python3 ../simple_modular_spacer_expansion_code/mapped_spacer_expansion_kmer_negative_control_in_spacer_in_parallel.py ../input_files/cas12a.fasta_all_hits.csv_genomes.fasta_crisprs.lst_reconciled_full_arr_positions.csv_no_blanks.csv ../input_files/cas12a.fasta_all_hits.csv_genomes.fasta_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_hitmap.csv_standardised.csv_non_redundant.csv_filtered.csv ../input_files/cas12a.fasta_all_hits.csv_genomes.fasta_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_hitmap.csv_standardised.csv_non_redundant.csv_filtered.csv_deduplicated.csv_filtered_hits_extracted_faidx_bp_window.fasta ../spacer_kmer_detection_all_systems/cas12a_neg_cntl_parallel_gadi_kmer14_29-3-2024.csv 14 52
