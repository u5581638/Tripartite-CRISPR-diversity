# script constructor 
# program to write the script required for GNU parallel

import sys
import csv
import os
#import file containing sequence names one on each line
#import file containing database file names one on each line

#extract the sequences
print(sys.argv)
block_directory = sys.argv[3]
perc_identity = sys.argv[4]
db_directory = sys.argv[5]
result_dir = sys.argv[6]
query_cover = sys.argv[7]
formatting = "\"10 qaccver saccver pident slen mismatch gapopen qcovs qlen sstart send evalue bitscore\""
#

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
if (os.path.isfile(db_directory + "parallelised_commands.sh")):
	os.remove(db_directory + "parallelised_commands.sh")
file = open (db_directory + "parallelised_commands.sh", "a")
for database in database_names:
	file.writelines([ "blastn ","-query ","\"" + db_directory + "spacer_distribution_analysis/"+ sys.argv[1].split("/")[-1] + "\"", " -db " + "\"" + block_directory + database.strip("\n") + "\""," -outfmt " + formatting,  " -perc_identity " + perc_identity,  " -max_target_seqs 100000000 -max_hsps 1 ",  "-qcov_hsp_perc " + query_cover + " -out \"" + db_directory + result_dir + sys.argv[1].split("/")[-1] + database.strip("\n") + "_all_hits.csv\"" + "\n" ])
# $db_directory$out_dir$sequence$database"_hits.csv"		
file.close()
	
