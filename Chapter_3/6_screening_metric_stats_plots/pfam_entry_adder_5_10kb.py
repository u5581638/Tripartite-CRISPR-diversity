# pfam_entry_adder

from Bio import SeqIO
import sys
import csv

# Note: This is a duplicate of pfam_entry_adder.py
# add Pfam scores to csv table containing co-occurrance scores for each representative protein sequence
def row_sorter(x):
	return float(x[3])

# INPUT: 
# 1. table of co-occurrance scores
# i.e. 0_10kb_full_co_occurrance_table_corrected_true_final_master.csv
with open(sys.argv[1],"r") as csvfile2:
	co_occurrance_table = list(csv.reader(csvfile2))

# 2. table containing Pfam predictions for each representative protein sequence (in csv format - layout generated by hhsuite_parser.py -> see Chapter 4).
# i.e. pfam_0_5kb_header_hits.csv
with open(sys.argv[2],"r") as csvfile:
	hit_table = list(csv.reader(csvfile))

# OUTPUT: table containing co-occurrance scores with the top pfam prediction appended.
# i.e. 0_10kb_full_co_occurrance_table_corrected_true_final_master.csv_pfam_added.csv
# SHELL: python3 pfam_entry_adder_5_10kb.py 0_10kb_full_co_occurrance_table_corrected_true_final_master.csv 

hit_dict = {}

for row in hit_table:
	if (row[0] not in hit_dict):
		hit_dict[row[0]] = [row]
	else:
		hit_dict[row[0]].append(row)
print(hit_dict.keys())
for row in hit_dict:
	hit_dict[row].sort(key=row_sorter,reverse=True)

ret_out = open(sys.argv[1] + "_pfam_added.csv","w")
spamwriter = csv.writer(ret_out)

for row in co_occurrance_table:
	if (row[0] in hit_dict):
		row.append(hit_dict[row[0]][0][1])
		row.append(hit_dict[row[0]][0][2])
		row.append(hit_dict[row[0]][0][3])
	spamwriter.writerow(row)

ret_out.close()

