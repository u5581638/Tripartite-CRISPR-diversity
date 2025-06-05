#!/bin/bash

#PBS -P do77
#PBS -q normalsr
#PBS -l walltime=48:00:00
#PBS -l mem=250GB
#PBS -l ncpus=52
#PBS -l wd
#PBS -l storage=scratch/xi88+gdata/xi88+gdata/va71

# is there any gpu specific flags

module load blast
module load python3/3.10.4
module load samtools

# perform spacer mapping and de-duplication against self-array mapping hits (for Type VI-D genomes)
# INPUT: 1. representative ORF/ or other CRISPR-associated protein sequence (in FASTA format)
# 		 i.e. cas13d.fasta (located in "queries/" folder)
# 		 2. file path to folder containing 10TB genome blocks, along with  files to run BLAST, for each block
#        i.e. /g/data/va71/labelled_genomes/
#		 3. file path to folder called "genomes/" containing the file with DNA extracted 20kb upstream and downstream of CRISPR-arrays
#		 i.e. "cas13d/genomes/" containing "spliced_debug_corrected_cleaned_rerun_20kb_local_windows_combined.fasta" and associated files for running BLAST
#		 4. Additionally, in order to run correctly, the output from running "" must be present in the "queries/" and "genomes/" directories
# OUTPUT: 1. Table containing hits to mapped spacers. This table has been filtered to exclude hits to self-CRISPR arrays, and hits to other arrays encodding the same spacers
#		  i.e. cas13d.fasta_all_hits.csv_genomes.fasta_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_hitmap.csv_standardised.csv
#		  2. FASTA file containing a list of "host" encoding contigs
#		  i.e. cas13d.fasta_all_hits.csv_genomes.fasta_crisprs.lst_full_real_arr_positions.csv
#		  3. FASTA file containing spacer mapped target contigs (in FASTA format)
#		  i.e. cas13d.fasta_all_hits.csv_genomes.fasta_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_hitmap.csv_standardised.csv_non_redundant.csv_filtered.csv_deduplicated.csv_filtered_hits_extracted_faidx_bp_window.fasta
#		  4. Preliminary table of PPS-spacer distances (not used)
#		  i.e. cas13d.fasta_all_hits.csv_genomes.fasta_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_hitmap.csv_standardised.csv_non_redundant.csv_filtered.csv_2_or_more_hits.csv_distances_annotated.csv
#		  5. Preliminary PAM concensous predictions for each subtype (not used)
#		  i.e. cas13d.fasta_all_hits.csv_genomes.fasta_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_hitmap.csv_standardised.csv_non_redundant.csv_filtered.csv_pam_summary.csv


echo $1
python3 /g/data/va71/crispr_pipeline_annotation/protein_annotation_tool_api_cli_in_progress_no_queries_mapping_draft_debug_full_spacer_control_mds_representatives_pfam_improved_SQL_spacer_only_v21.py -i /g/data/va71/crispr_pipeline_annotation/annotation_upgraded_main_workflow_run/cas13d/queries/cas13d.fasta -d /g/data/va71/crispr_pipeline_annotation/annotation_upgraded_main_workflow_run/cas13d/genomes/ -t "query" -m "/g/data/va71/labelled_genomes/" -b -r -c -tony -orientation -crispr_detect -phage_genomes -PAM_anno -phage_rep -spacer_gen_bypass -pps_map -prot_case_generation -cores 52 &> spacer_mapping_cas13b_dump.txt