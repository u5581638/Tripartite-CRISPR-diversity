#!/bin/bash

#PBS -P do77
#PBS -q normalsr
#PBS -l walltime=48:00:00
#PBS -l mem=52GB
#PBS -l ncpus=26
#PBS -l wd
#PBS -l storage=scratch/xi88+gdata/xi88+gdata/va71

# is there any gpu specific flags

module load blast
module load python3/3.10.4
module load samtools

# INPUT: representative ORF/ or other CRISPR-associated protein sequence (in FASTA format)
# also requires the folder the 10TB database genome blocks are located in. This is vestigal.
# also requires the files containing the 40kb windows (DNA extracted up to 20kb upstream and downstream of the CRISPR-array)
# OUTPUT: Set of contigs containing the query sequence within 20kb of a CRISPR array. A table describing a non-redundant set of predicted CRISPR-arrays by 3 seperate CRISPR-array prediction tools. Overlapping predictions were ranked (CRISPRdetect > CRT > PILERCR) by their overall accurary

echo $1
python3 /g/data/va71/crispr_pipeline_annotation/protein_annotation_tool_api_cli_in_progress_no_queries_mapping_draft_debug_full_spacer_control_mds_representatives_pfam_improved_SQL_premapping_v21.py -i /g/data/va71/crispr_pipeline_annotation/annotation_upgraded_main_workflow_run/cas13a/queries/cas13a.fasta -d /g/data/va71/crispr_pipeline_annotation/annotation_upgraded_main_workflow_run/cas13a/genomes/ -t "query" -m "/g/data/va71/labelled_genomes/" -b -r -c -tony -orientation -rnafold_switch -crispr_detect -phage_genomes -PAM_anno -phage_rep -prot_case_generation -cores 26 &> prespacer_cas13a_dump.txt