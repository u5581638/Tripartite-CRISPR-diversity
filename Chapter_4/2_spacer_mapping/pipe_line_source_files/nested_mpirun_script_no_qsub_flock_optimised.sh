#!/bin/bash

#PBS -P do77
#PBS -q normalsr
#PBS -l walltime=48:00:00
#PBS -l mem=24GB
#PBS -l ncpus=12
#PBS -l wd
#PBS -l storage=scratch/xi88+gdata/xi88+gdata/va71


module load blast/2.10.1
module load parallel
module load python3/3.10.4


# nested GNU parallel script
# script to call GNU parallel for each sequence which in turn individually calls GNU parallel for database block.
# results should be manually concatenated together at the end of each run to avoid truncation.
# top level script need only be an GNU parallel if required/ to increase speed
# easily convertable to a self-submitting job script.

# each sequence script should:
# 1. generate an individual script for running GNU parallel
# 2. run GNU parallel from a given directory
# 3. concatenate the results

# input: script containing the input sequence file names listed one per line
# these should be readable as csv

spacer_input=$1
block_directory=$2
perc_identity=$3
cores=$4
working_dir=$5
db_dir=$6
arr=$spacer_input
spacer_result_dir=$7
query_cover=$8
formatting=$9


export ncores_per_task=1

 

if [ -f "/g/data/va71/crispr_pipeline_annotation/pipe_line_source_files/nucleotide_database_headers.txt" ] 
then 
	python3 /g/data/va71/crispr_pipeline_annotation/pipe_line_source_files/script_constructor_all_calls_direct.py $spacer_input /g/data/va71/crispr_pipeline_annotation/pipe_line_source_files/nucleotide_database_headers.txt $block_directory $perc_identity $db_dir $spacer_result_dir $query_cover $formatting  #creates the file: single_parallelised_commands.sh
else 
	find $block_directory -name "*.fasta" -type "f" -exec basename {} \; |  xargs -n 1 -I {} -P 1 cat {} >> nucleotide_database_headers.txt # check this line of code on Monday to see if it works
	python3 /g/data/va71/crispr_pipeline_annotation/pipe_line_source_files/script_constructor_all_calls_direct.py $spacer_input /g/data/va71/crispr_pipeline_annotation/pipe_line_source_files/nucleotide_database_headers.txt $block_directory $perc_identity $db_dir $spacer_result_dir $query_cover $formatting
fi 


touch ${spacer_input}"_all_hits.csv"
flock ${spacer_input}"_all_hits.csv" -c  'parallel -P $cores -a ${db_dir}parallelised_commands.sh'
 
# touch "../run_complete.txt"
find ${db_dir}${spacer_result_dir} -name "*_hits.csv" -type "f" | xargs -n 1 -I {} -P 1 cat {} >> ${spacer_input}"_all_hits.csv"

rm ${db_dir}parallelised_commands.sh
# rm -r ${db_dir}${spacer_result_dir}

# the previous 4 lines could easily be placed in a single script and qsubbed alone.
