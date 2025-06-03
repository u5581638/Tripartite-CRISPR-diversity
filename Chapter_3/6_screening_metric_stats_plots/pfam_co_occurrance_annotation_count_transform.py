# script to transform the representation of the co_occurrance annotations to make counting faster
# as well as labelling the points.

import csv
import sys


def row_sorter(x):
	return x[2]

# INPUT: Table of representative protein sequences with annotation classification appended.
# i.e. 0_10kb_full_co_occurrance_table_corrected_true_final_master.csv_pfam_added.csv
# OUTPUT:
# i.e. 0_10kb_full_co_occurrance_table_corrected_true_final_master.csv_padloc_added.csv_pfam_count.csv
# SHELL: python3 padlocplus_co_occurrance_annotation_count_transform.py 0_10kb_full_co_occurrance_table_corrected_true_final_master.csv_pfam_added.csv

with open(sys.argv[1],"r") as csvfile:
	hit_table = list(csv.reader(csvfile))

# want to be able to label genomes.
# build a dictionary with the annotation tool as the key.
# values should be [annotation_id,count,[list of sequence identifers]]
# then write values to table.

hit_dict = {}
for row in hit_table[1:]:
	if (row[10] not in hit_dict):
		hit_dict[row[10]] = [row[10],row[11], 1, [row[0]]]
	else:
		hit_dict[row[10]][2] += 1
		hit_dict[row[10]][3].append(row[0])

out_results = list(hit_dict.values())
out_results.sort(key=row_sorter,reverse=True)
ret_out = open(sys.argv[1] + "_pfam_count.csv","w")
spamwriter = csv.writer(ret_out)

for row in out_results:
	spamwriter.writerow(row)

ret_out.close()	