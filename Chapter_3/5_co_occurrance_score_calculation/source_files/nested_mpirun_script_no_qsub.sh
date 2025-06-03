#!/bin/bash

#PBS -P do77
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -l mem=96GB
#PBS -l ncpus=48
#PBS -l wd
#PBS -l storage=scratch/xi88+gdata/xi88+gdata/va71


module load blast/2.10.1
module load parallel
module load python3/3.8.5
pip3 install biopython

# nested tblastn script
# script to call tblastn for each sequence which in turn individually calls tblastn for database block.
# results should be manually concatenated together at the end of each run to avoid truncation.
# top level script need only be an tblastn if required/ to increase speed
# easily convertable to a self-submitting job script.

# each sequence script should:
# 1. generate an individual script for running tblastn
# 2. run tblastn from a given directory
# 3. concatenate the results

# INPUT: 1. folder containing containing the input sequences named as run*.fasta
# 		 2. Folder must also contain nucleotide_database_headers.txt
#		 3. although not a direct input, each genome block in labelled/genomes must also be available to perform tBLASTn
# OUTPUT: list of tBLASTn hits from each genome block in co_occurrance_bug_diagnostics_results.  These may be concatenated into one result file.
# i.e. combined_run_1.csv
# these should be readable as csv

arr=($( find ../ -name "run_*.fasta" -type "f" -exec basename {} \; )) # need to change directory to a new folder.

export ncores_per_task=1
for ele in ${arr[*]} 
do 
echo $ele >> sequence_header.txt
done

python3 script_constructor_all_calls.py sequence_header.txt nucleotide_database_headers.txt #creates the file: single_parallelised_commands.sh
mkdir ../co_occurrance_bug_diagnostics_results/
parallel -P 48 -a parallelised_commands.sh
mkdir ../co_occurrance_bug_diagnostics_final_results 

#find ../co_occurrance_bug_diagnostics_results/ -name "*_hits.csv" -type "f" | xargs -n 1 -I {} -P 1 >> "../co_occurrance_bug_diagnostics_final_results/combined_"${arr[0]}".csv"

touch "../run_complete.txt"
rm -r ../touch_files
mkdir ../touch_files

# the previous 4 lines could easily be placed in a single script and qsubbed alone.
