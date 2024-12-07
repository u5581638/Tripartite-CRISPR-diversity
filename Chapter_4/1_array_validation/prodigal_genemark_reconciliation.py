# This file contains two functions to take either the union or intersection of prodigal and Genemark predicted genes.

import sys
from Bio import SeqIO

def dash_inserter (my_list):
	ret_str = my_list[0]
	i = 1
	while (i < len(my_list)):
		ret_str = ret_str + "_" + my_list[i]
		i += 1
	return ret_str	

def union (genemark_orfs_url, prodigal_orfs_url, protein_out_url): # take every non-redundant orf

	genemark_sequences = SeqIO.parse(genemark_orfs_url, "fasta")
	prodigal_sequences = SeqIO.parse(prodigal_orfs_url, "fasta")
	seq_dict = {}
	for sequence in genemark_sequences:
		genome_id = dash_inserter(sequence.id.split("_") [:-1] )
		genome_description = sequence.description.split(" # ")
		protein_start = genome_description[1]
		protein_end = genome_description[2]
		protein_sense = genome_description[3]
		if (protein_sense == "+"):
			protein_sense = "1"
		elif (protein_sense == "-"):
			protein_sense = "-1"
		else:
			print("error")
			print(protein_sense)
			exit()		
		if ((protein_start, protein_end, protein_sense) not in seq_dict):
			seq_dict[(protein_start, protein_end, protein_sense)] = sequence

	for sequence in prodigal_sequences:
		genome_id = dash_inserter(sequence.id.split("_") [:-1] )
		genome_description = sequence.description.split(" # ")
		protein_start = genome_description[1]
		protein_end = genome_description[2]
		protein_sense = genome_description[3]
		if ((protein_start, protein_end, protein_sense) not in seq_dict):
			seq_dict[(protein_start, protein_end, protein_sense)] = sequence

	out_values = seq_dict.values()
	SeqIO.write(out_values, protein_out_url, "fasta")
	return out_values

def intersection (genemark_orfs_url, prodigal_orfs_url, protein_out_url):
	genemark_sequences = SeqIO.parse(genemark_orfs_url, "fasta")
	prodigal_sequences = SeqIO.parse(prodigal_orfs_url, "fasta")
	genemark_dict = {}
	for sequence in genemark_sequences:
		genome_id = dash_inserter(sequence.id.split("_") [:-1] )
		genome_description = sequence.description.split(" # ")
		protein_start = genome_description[1]
		protein_end = genome_description[2]
		protein_sense = genome_description[3]
		if (protein_sense == "+"):
			protein_sense = "1"
		elif (protein_sense == "-"):
			protein_sense = "-1"
		else:
			print("error")
			print(protein_sense)
			exit()		
		genemark_dict[(protein_start, protein_end, protein_sense)] = sequence
	prodigal_dict = {}
	for sequence in prodigal_sequences:
		genome_id = dash_inserter(sequence.id.split("_") [:-1] )
		genome_description = sequence.description.split(" # ")
		protein_start = genome_description[1]
		protein_end = genome_description[2]
		protein_sense = genome_description[3]
		prodigal_dict[(protein_start, protein_end, protein_sense)] = sequence
	ret_list = []
	for my_key in prodigal_dict:
		if my_key in genemark_dict:
			ret_list.append(genemark_dict[my_key])	
	SeqIO.write(ret_list, protein_out_url, "fasta")		
	return ret_list
