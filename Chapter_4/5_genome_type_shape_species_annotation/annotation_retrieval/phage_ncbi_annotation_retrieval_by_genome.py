import csv
import sys


# script to retrieve all the ncbi annotation for a given cluster input file:
# should have an existing GOLD equivalent: "phage_gold_annotation_retrieval.py"
with open(sys.argv[1],"r") as csvfile:
	cluster_table = list(csv.reader(csvfile))

all_matching_gold_anno = []

for cluster in cluster_table[1:]:
	my_entries = [cluster[0]]
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

# outgold_url = open(sys.argv[2],"w")
outgold_url = open(sys.argv[2],"w")
spamwriter = csv.writer(outgold_url)

# This loop could be parallalised!!
for entry in block_id_dict.keys():
	# may need to replace with sys.argv[2]
	print("entry:")
	print(entry) 
	seq_url = "/g/data/va71/labelled_genomes/" + entry + ".fasta_gold_annotations.csv_missed.csv_ncbi.csv" # need to fill in the exact url
	seq_url2 = "/g/data/va71/labelled_genomes/" + entry + ".fasta_all_ncbi.csv"
#	seq_url = sys.argv[2]
#	seq_url2 = sys.argv[3]
	# need to load both files
	with open(seq_url,"r") as csvfile:
		block = list(csv.reader(csvfile)) # may be able to not load this in memory
	
# may need to throw an exception if this fails!!
	with open(seq_url2,"r") as csvfile:
		block2 = list(csv.reader(csvfile)) # may be able to not load this in memory
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
	for sequ in block2:
		print("second_block:")
		print(entry + "|" + sequ[0])
		if (entry + "|" + sequ[0] in my_id_set):
			spamwriter.writerow([entry + "|" + sequ[0],"","","","","","","",sequ[1]])
			pos += 1
		else:
			neg += 1
				

#	print("Proportion2:")
	


outgold_url.close()