# ncbi_window_microbiome_deduplication

#check contig deduplcation has not inflated ncbi taxonomy estimations

import csv
import sys

with open(sys.argv[1],"r") as csvfile:
	ncbi_table = list(csv.reader(csvfile))

ncbi_dict = {}
for row in ncbi_table:
	if (row[3] not in ncbi_dict):
		ncbi_dict[row[3]] = row
ret_out = open(sys.argv[2],"w")
spamwriter = csv.writer(ret_out)
for row in ncbi_dict:
	spamwriter.writerow(ncbi_dict[row])
ret_out.close()
