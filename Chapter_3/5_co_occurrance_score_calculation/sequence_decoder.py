# program to split the results file containing many blast hits into individual arbitarily named results files.
import csv 
import sys

table_dict = {} # dictionary containing each sequence_id and a list of 

# INPUT: csv file containing search results from tBLASTn search using multiple queries of representative protein sequences derived from clustering
# i.e. combined_run_1.fasta_all_hits.csv
# i.e  spliced_cleaned_tblastn_windows.fasta_local_window_hits.csv (need actual file name)
# OUTPUT: list of aggregated tBLASTn results for each protein sequence (at either local window or global level, depending on the input file)
# i.e. sequence_tmp_sequence_1_combined_run_1.fasta_all_hits.csv
# SHELL: python3. sequence_decoder.py combined_run_1.fasta_all_hits.csv

with open (sys.argv[1], "r") as csvfile:
	hit_table = csv.reader(csvfile)

	row_num = 0
	for row in hit_table:
		if (row[0] not in table_dict):

			print("new_row")
			print(row_num)
			row_num += 1
			table_dict[row[0]] = [row]
		else: # row in hit table
			table_dict[row[0]].append(row)	

i = 0
print("Ready to write!!")
for sequence_results in table_dict:
	print(sequence_results)
	sequence_id = sequence_results
	sequence_id = sequence_id.split("|")
	sequence_id = sequence_id[0] + sequence_id[1]
#	sequence_id = sequence_id.replace("|","_")
	print(sequence_id)
	out_file = open("sequence_tmp_" + sequence_id + "_" + sys.argv[1], "a")
	spam_writer = csv.writer(out_file)
	for row in table_dict[sequence_results]:
		spam_writer.writerow(row)
	i += 1
	out_file.close()	




