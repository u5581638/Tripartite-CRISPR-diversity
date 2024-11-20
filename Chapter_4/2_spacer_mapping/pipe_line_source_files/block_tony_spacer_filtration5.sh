#!/bin/bash
# -P va71
# -q normal
# -l walltime=48:00:00
# -l mem=190GB
# -l ncpus=6
# -l wd
# -l storage=scratch/xi88+gdata/xi88+gdata/va71

# batch job to filter spacers from crispr-self-targetted hits by block
module load python3/3.10.4

working_dir=$1
mkdir ${working_dir}"co_occurrance_bug_diagnostics_filtered_results/"
find  ${working_dir}co_occurrance_bug_diagnostics_results2/ -name "*" -type "f" -exec basename {} \; | xargs -n 1 -I {} -P 6 /g/data/va71/crispr_pipeline_annotation/pipe_line_source_files/block_tony_filtration_runner.sh {} $working_dir 5