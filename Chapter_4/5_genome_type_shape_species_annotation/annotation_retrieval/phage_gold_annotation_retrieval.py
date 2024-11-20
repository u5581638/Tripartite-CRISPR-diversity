# phage_gold_sequence_retrieval

# script to retrieve gold annotations for spacer mapped phage sequences

# input = cluster file of mapped phage sequences
# 		= labelled genomes gold files

# Approach:
# 1. Sort the identifers in alphabetical order.
# 2. Load genome once per each set of identifiers
# 3. Write matching annotations to output file.


import csv
import sys


with open(sys.argv[1],"r") as csvfile:
	cluster_table = list(csv.reader(csvfile))

all_matching_gold_anno = []



for cluster in cluster_table[1:]:
	my_entries = cluster[7].split(" ")
	all_matching_gold_anno.extend(my_entries)

block_id_dict = {}

for entry in all_matching_gold_anno:
	my_entry = entry.split("|")
	block_id = my_entry[0]
	full_entry = entry.split(":") [0]
	if (block_id not in block_id_dict):
		block_id_dict[block_id] = [full_entry]
	else:
		block_id_dict[block_id].append(full_entry)			

outgold_url = open(sys.argv[2],"w")
spamwriter = csv.writer(outgold_url)

# This loop could be parallalised!!
for entry in block_id_dict.keys():
	seq_url = "/g/data/va71/labelled_genomes/" + entry + ".fasta_gold_annotations.csv" # need to fill in the exact url
	with open(seq_url,"r") as csvfile:
		block = list(csv.reader(csvfile)) # may be able to not load this in memory
	my_id_set = set(block_id_dict[entry])
	pos = 0
	neg = 0
	# the outfile may be reduced by taking only the set of matches.
	for sequ in block:
		if (sequ[0] in my_id_set):
			spamwriter.writerow(sequ)
			pos += 1
		else:
			neg += 1
	# removed this as the proprotion is meaningless
	# print("Proportion:")
	# print((pos / neg))	


outgold_url.close()




