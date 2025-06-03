# ncbi_window_microbiome_deduplication

#check contig deduplcation has not inflated ncbi taxonomy estimations

import csv
import sys

# INPUT: ncbi table of species annotations in csv format.
# i.e. ncbi_microbiome_window_only2.csv
# OUTPUT: jgi/ncbi table with entries from the same contig deduplicated
# i.e. ncbi_microbiome_window_only2.csv_dedup.csv
# SHELL: ncbi_window_microbiome_deduplication.py ncbi_microbiome_window_only2.csv ncbi_microbiome_window_only2.csv_dedup.csv

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
