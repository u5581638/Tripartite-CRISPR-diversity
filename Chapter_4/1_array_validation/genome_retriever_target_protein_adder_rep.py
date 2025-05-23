# genome_retriever_target_protein_adder


import sys
from Bio import SeqIO
import csv
import os

# retrieve the corresponding genomes and annotate them with the target protein sense and position coordinates
def protein_adder (sequence_url, genomes_st):
	# prodigal predicted genomes
	sequences = SeqIO.parse(sequence_url, "fasta")
	# represenative geonomes
	genomes_handle = genomes_st
	genomes = {}
	for genome in genomes_handle:
		genomes[genome.id] = genome


	ret_genomes = []
	for sequ in sequences:
		my_sequ_id = sequ.id
		my_sequ_id_genome_split_coord = my_sequ_id.split("::")
		my_sequ_base = my_sequ_id_genome_split_coord[0]
		my_sequ_id_genome_split_coord = my_sequ_id_genome_split_coord[1]
		my_sequ_id_genome_split_coord = my_sequ_id_genome_split_coord.split("_")
		my_sequ_order = my_sequ_id_genome_split_coord[1]
		my_sequ_id_genome_split_coords = my_sequ_id_genome_split_coord[0]
		complete_id = my_sequ_base + "::" + str(my_sequ_id_genome_split_coords)

		my_sequ_params = sequ.description.split(" # ")
		my_sequ_target_start = my_sequ_params[1]
		my_sequ_target_end = my_sequ_params[2]
		my_sequ_target_sense = my_sequ_params[3]

		if (complete_id in genomes):
			my_genome = genomes[complete_id]
			my_genome.description = my_genome.description + " " + "tar_prot_start_pos-" + my_sequ_target_start + " " + "tar_prot_end_pos-" + my_sequ_target_end + " " + "tar_prot_sense-" + my_sequ_target_sense +" " + "tar_prot_order-" + my_sequ_order
			ret_genomes.append(my_genome)
		else:
			print("Bugger!!")
			print(complete_id)	
	return ret_genomes	