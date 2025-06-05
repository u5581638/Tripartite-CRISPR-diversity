#!/bin/bash

#PBS -P do77
#PBS -q normalsr
#PBS -l walltime=48:00:00
#PBS -l mem=250GB
#PBS -l ncpus=52
#PBS -l wd
#PBS -l storage=scratch/xi88+gdata/xi88+gdata/va71

module load python3/3.10.4
module load samtools
module load blast

# input urls to run virsorter2 only for the mapped sequences of type VI-A systems
# INPUT: 1. path to folder containing DNA extracted within 20kb of CRISPRs (db_directory_path)
#		 i.e. /g/data/va71/crispr_pipeline_annotation/annotation_upgraded_main_workflow_run/cas13a/genomes/
#		 2. query search name 
#		 i.e. cas13a.fasta
#		 3. tBLASTn result table from contig retrieval.
#		 i.e. /g/data/va71/crispr_pipeline_annotation/annotation_upgraded_main_workflow_run/cas13a/queries/cas13a.fasta_all_hits.csv

# OUTPUT: folder containing virsorter predictions for host/mapped contigs

db_dir=/g/data/va71/crispr_pipeline_annotation/annotation_upgraded_main_workflow_run/cas13a/genomes/
b_basename=cas13a.fasta
all_hits=/g/data/va71/crispr_pipeline_annotation/annotation_upgraded_main_workflow_run/cas13a/queries/cas13a.fasta_all_hits.csv
cores=52
virsort=1

python3 /g/data/va71/crispr_pipeline_annotation/annotation_parallelisation_standalone_running_cmds_plasme_general_all.py $db_dir $b_basename $all_hits $cores $virsort
