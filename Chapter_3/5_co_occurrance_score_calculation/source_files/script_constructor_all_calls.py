# script constructor 
# program to write the script required for nci-parallel

import sys
import csv

#import file containing sequence names one on each line
#import file containing database file names one on each line

#extract the sequences
sequence_names = []
with open(sys.argv[1], "r") as seq_file:
	while (True):
		line = seq_file.readline()
		if ("" == line):
			break
		sequence_names.append(str(line))

#extract the database names
database_names = []
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
	
