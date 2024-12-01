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

# input urls to run virsorter2 only for the mapped sequences of type VI-D systems


db_dir=/g/data/va71/crispr_pipeline_annotation/annotation_upgraded_main_workflow_run/cas13d/genomes/
b_basename=cas13d.fasta
all_hits=/g/data/va71/crispr_pipeline_annotation/annotation_upgraded_main_workflow_run/cas13d/queries/cas13d.fasta_all_hits.csv
cores=52
virsort=1

python3 /g/data/va71/crispr_pipeline_annotation/annotation_parallelisation_standalone_running_cmds_plasme_general_all.py $db_dir $b_basename $all_hits $cores $virsort
