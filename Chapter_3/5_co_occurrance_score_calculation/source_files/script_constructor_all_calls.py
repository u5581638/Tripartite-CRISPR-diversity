# script constructor 
# helper program to write the script required for parallel tblastn search

import sys
import csv

# INPUT: 1. import file containing protein sequence file names one on each line
#        i.e. sequence_headers.txt
# 		 2. import file containing database block file names one on each line
#		 i.e. nucleotide_database_headers.txt
# OUTPUT: BASH file containing TBLASTN commmands for each genome block.
# SHELL: N/A (called from nested_mpirun_script_no_qsub.sh)
# extract the sequences
sequence_names = []

with open(sys.argv[1], "r") as seq_file:
	while (True):
		line = seq_file.readline()
		if ("" == line):
			break
		sequence_names.append(str(line))

# extract the database names
database_names = []
# import file containing database block file names one on each line
with open(sys.argv[2], "r") as database_names_f:
	while (True):
		line = str(database_names_f.readline())
		if ("" == line):
			break
		database_names.append(line)

file = open ("parallelised_commands.sh", "a")
for sequ in sequence_names:
	
	for database in database_names:
		file.writelines(["./parallel_pipeline_runner.sh "," ../" + sequ.strip("\n") , " /g/data/va71/labelled_genomes/" + database.strip("\n") + "\n" ])
file.close()
	
