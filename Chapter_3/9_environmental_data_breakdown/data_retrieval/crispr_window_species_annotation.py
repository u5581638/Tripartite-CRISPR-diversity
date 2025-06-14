# crispr_window_species_annotation

from Bio import SeqIO
import sys
import csv

# Note: Somewhat misleadingly, this file actually finds gold annotations for the 10TB data block, rather than the individual windows.
# script to retrieve sequence metadata for JGI-derived sequences using a merged metadata table downloaded from the GOLD database website

# SHELL: find labelled_genomes -name "*.fasta" | xargs -n 1 -I {} -P 1 python3 crispr_window_species_annotation.py {} true_gold_compile3_merged.csv {}_gold_annotations.csv
# INPUT: 1. genome "blocks" in FASTA format.
#		 2. merged_table containing JGI sequence identifers (Ga numbers). This must be: true_gold_compile3_merged.csv. Any other table will not work!!
# OUTPUT: table containing matching GOLD annotations.

sequences = SeqIO.parse(sys.argv[1],"fasta")

with open(sys.argv[2], "r") as csvfile:
	hit_table = csv.reader(csvfile)
	hit_dict = {}
	for row in hit_table:
		if (row[0] not in hit_dict):
			hit_dict[row[0]] = row
		else:
			print("Error!!")


ret_out = open(sys.argv[3],"w")
spamwriter = csv.writer(ret_out)
spamwriter.writerow(["Sequence_writer","AP_id","AP_NAME","AP_GENBANK","AP_SRA","AP_GOLD_PROJECT_ID","PROJECT_NAME","SEQUENCING_STRATEGY","STUDY_ID","ORGANISM_ID","BIOSAMPLE_ID","ORGANISM_NAME","PHYLUM","CLASS","ORDER","FAMILY","GENUS","SPECIES","GRAM","BIOSAMPLE_NAME","NAME","COLLECTION_SITE","LOCATION","BIOSAMPLE_ECOSYSTEM","ECOSYSTEM_TYPE","ECOSYSTEM_SUBTYPE","ECOSYSTEM"])
missed_list = open(sys.argv[3] + "_missed.csv","w")
missed_writer = csv.writer(missed_list)
for sequ in sequences:
	sequ_id = sequ.id.split("|") [1]
	sequ_id = sequ_id.split("_") [0]
#	print(sequ_id)
#	print(hit_dict.keys())
#	exit()
	if sequ_id in hit_dict:
		spamwriter.writerow([sequ.id] + hit_dict[sequ_id])
	else:
		missed_writer.writerow([sequ.id])




ret_out.close()
missed_list.close()