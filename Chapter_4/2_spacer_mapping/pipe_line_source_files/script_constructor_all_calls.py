# script constructor 
# program to write the script required for GNU parallel

import sys
import csv

#import file containing sequence names one on each line
#import file containing database file names one on each line

#extract the sequences
print(sys.argv)
block_directory = sys.argv[3]
perc_identity = sys.argv[4]
db_directory = sys.argv[5]
result_dir = sys.argv[6]
query_cover = sys.argv[7]
formatting = sys.argv[8]


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

# check this if this file exists first!!!
file = open (db_directory + "parallelised_commands.sh", "a")
for sequ in sequence_names:
	for database in database_names:

		# 
		
		file.writelines(["/g/data/va71/crispr_pipeline_annotation/pipe_line_source_files/" + "parallel_pipeline_runner.sh ",sequ.strip("\n") , " " + block_directory + database.strip("\n")," " + perc_identity,  " " + db_directory,  " " + result_dir, " " + query_cover + " " + formatting + " " + "\n" ])
file.close()
	
