# table_identifier_extractor

import csv
import sys
from Bio import SeqIO
import copy

# function to filter hits (from a BLAST search with output in csv format) based on the presence of matching genome sequence identifiers (from a file containing genomes in FASTA format). These filtered hits are then written to a seperate file
def hits_from_genomes (genomes, table_url, r=False):
	
	with open(table_url, "r") as csvfile:
		hit_table = list(csv.reader(csvfile))	
	gene_dict = {}
	for gene in genomes:
		my_description = copy.deepcopy(gene.description)
		my_description_id = gene.id
		my_description = my_description.split("|")
		my_description = my_description[1]
		if (r == True):
			my_description = my_description.split(" ")
			my_description = my_description[0]
			my_id = my_description + "|" + my_description_id
		else:
			my_id = my_description_id	

		gene_dict[my_id] = gene
	ret_url = open(table_url + "_rep_genome.fasta_hits.csv", "w")
	spam_writer = csv.writer(ret_url)
	for hit in hit_table:
	#	print(hit)
		hit_id = hit[1].split("|")
		hit_id = hit_id[1]
	#	print(hit_id)
		if (hit_id in gene_dict):
			spam_writer.writerow(hit)
			print("Target_found!!")
		else:
			pass

	ret_url.close()
	return table_url + "_rep_genome.fasta_hits.csv"				



