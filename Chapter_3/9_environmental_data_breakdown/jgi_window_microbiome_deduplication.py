# ncbi_window_microbiome_deduplication

# ncbi taxonomy estimations without contig duplication

import csv
import sys

# INPUT: jgi/ncbi table of species annotations in csv format.
# i.e. jgI_microbiome_window_only2.csv
# OUTPUT: jgi/ncbi table with entries from the same contig deduplicated
# i.e. jgI_microbiome_window_only2.csv_dedup.csv
# SHELL: jgi_window_microbiome_deduplication.py jgI_microbiome_window_only2.csv jgI_microbiome_window_only2.csv_dedup.csv

with open(sys.argv[1],"r") as csvfile:
	ncbi_table = list(csv.reader(csvfile))

ncbi_dict = {}
for row in ncbi_table:
	if (row[12] not in ncbi_dict):
		ncbi_dict[row[12]] = row

ret_out = open(sys.argv[2],"w")
spamwriter = csv.writer(ret_out)
for row in ncbi_dict:
	spamwriter.writerow(ncbi_dict[row])
ret_out.close()
