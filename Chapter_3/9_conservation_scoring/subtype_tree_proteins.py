# program to extract sequences from a given CRISPR-subtype/all subtypes for tree generation with a given PFAM annotation.

import sys
import csv
from Bio import SeqIO

# CRISPR_subtype annotation_table
csvfile = open(sys.argv[1],"r")
annotation_table = csv.reader(csvfile)

sequences = SeqIO.parse(sys.argv[2],"fasta")

seq_dict = {}

for sequ in sequences:
	if sequ.id not in seq_dict:
		seq_dict[sequ.id] = sequ 

query = sys.argv[3]
query_index = 14 # column containing pfam/padloc entries
# specify output
ret_seqs = open ( sys.argv[1] + "_" + query + "_phylo_seqs.fa","a")
spamwriter = csv.writer(ret_seqs)
for row in annotation_table:
	profile_id = row[query_index]
	my_id = profile_id.split(" ") [0]
	if my_id == query:
		protein_index_num = row[1].split("_")[-1]
		protein_id = row[1].split("::")
		protein_range = protein_id[1].split("_")[0]
		profile_id = protein_id[0]
		protein_index = profile_id + "::" + protein_range + "_" + protein_index_num
		if (protein_index in seq_dict):
			SeqIO.write(seq_dict[protein_index],ret_seqs,"fasta")
		else:
			print("Error, could not find protein:",protein_index)

csvfile.close()
ret_seqs.close()
