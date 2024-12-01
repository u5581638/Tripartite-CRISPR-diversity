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

# perform spacer mapping and de-duplication against self-array mapping hits (for Type VI-B genomes)

echo $1
python3 /g/data/va71/crispr_pipeline_annotation/protein_annotation_tool_api_cli_in_progress_no_queries_mapping_draft_debug_full_spacer_control_mds_representatives_pfam_improved_SQL_spacer_only_v21.py -i /g/data/va71/crispr_pipeline_annotation/annotation_upgraded_main_workflow_run/cas13b/queries/cas13b.fasta -d /g/data/va71/crispr_pipeline_annotation/annotation_upgraded_main_workflow_run/cas13b/genomes/ -t "query" -m "/g/data/va71/labelled_genomes/" -b -r -c -tony -orientation -crispr_detect -phage_genomes -PAM_anno -phage_rep -spacer_gen_bypass -pps_map -prot_case_generation -cores 52 &> spacer_mapping_cas13b_dump.txt