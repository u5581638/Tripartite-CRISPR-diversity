# co_occurance_calculator
# program to calculate the co_occurrance score for a given protein given the sequence

import sys
import csv

# INPUT:
# 1. argv[1] = protein sequence
# 2. argv[2] = csv file containing number of hits inside CRISPR array
# 3. argv[3] = csv file containing of database hits (with \n and headers removed)
# 4. argv[4] = result table tabulating the co-occurrance scores

# OUTPUT: co-occurrance score for the protein sequence used as input.
# SHELL: This program is usually run as a helper function for co_occurrance_batch_helper_script.sh, which is in turn called by co_occurrance_calculation_script.sh

#print (sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
with open(sys.argv[2], "r") as csvfile:
	window_hit_table = list(csv.reader(csvfile))
with open(sys.argv[3], "r") as csvfile2:
	database_hit_table = list(csv.reader(csvfile2))
	#print("program started")
with open(sys.argv[4], "a") as csvfile3:
	spam_writer = csv.writer(csvfile3)
	spam_writer.writerow([sys.argv[1], sys.argv[2], sys.argv[3], len(window_hit_table), len(database_hit_table), float(len(window_hit_table) / len(database_hit_table))])
	
