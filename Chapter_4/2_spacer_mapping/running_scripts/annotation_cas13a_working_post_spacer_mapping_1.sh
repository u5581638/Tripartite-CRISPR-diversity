#!/bin/bash

#PBS -P do77
#PBS -q normalsr
#PBS -l walltime=48:00:00
#PBS -l mem=96GB
#PBS -l ncpus=12
#PBS -l wd
#PBS -l storage=scratch/xi88+gdata/xi88+gdata/va71

# is there any gpu specific flags

module load blast
module load python3/3.10.4
module load samtools

# generate PPS-distance and other files after spacer mapping yet prior to analysing the gene composition (for Type VI-A genomes)
# uses the same inputs as the original spacer mapping script
# In cases where there is not sufficent compute to run as a single continuous batch job
echo $1
python3 /g/data/va71/crispr_pipeline_annotation/protein_annotation_tool_api_cli_in_progress_no_queries_mapping_draft_debug_full_spacer_control_mds_representatives_pfam_improved_SQL_spacer_only_v21.py -i /g/data/va71/crispr_pipeline_annotation/annotation_upgraded_main_workflow_run/cas13a/queries/cas13a.fasta -d /g/data/va71/crispr_pipeline_annotation/annotation_upgraded_main_workflow_run/cas13a/genomes/   -t "query" -m "/g/data/va71/labelled_genomes/" -b -r -c -orientation -PAM_anno -pps_map -crispr_detect -phage_genomes -spacer_gen_bypass -ignore_hitmap_sql -filter_hits -prot_case_generation -ignore_initialisation -rnafold_switch -skip -cores 12 &> stdout_sterr_cas13a_post_mapping_dump1.txt