# CRISPR array position appender

# the plan.

# load relabelled file
# load concatenated CRISPR-array predictions produced from PILER-CR
# extract position score
# Find the minimum distance between the array and the starting point of the arrays
# this requires computing the abolsute positions of the proteins from the CRISPR array
# the abolsute position will be the position given by prodigal/genemark + starting array position
# if strand is negative the starting position should have a higer score than the negative position
# multiple mapping between arrays and ids mean that the minimal distance between the ids and arrays must be found

from Bio import SeqIO
import csv
import sys
import re

# filtered cluster files of protein sequences in FASTA format.
sequences = SeqIO.parse(sys.argv[1], "fasta")
# CRISPR-array predictions from PILER-CR in FASTA format.
crispr_arrs = list(SeqIO.parse(sys.argv[2], "fasta"))


crispr_arr_dict = {}

def closest_crispr_array (arr_positions, protein_position, first_gene_coord):
	
	true_position = protein_position + first_gene_coord # true location of protein
	nearest_position= abs(arr_positions[0] - true_position)
	nearest_arr_position = arr_positions[0]
	for position in arr_positions:
		if (abs(position - true_position) <  nearest_position):
			nearest_position = abs(position - true_position)
			nearest_arr_position = position	
	return (nearest_arr_position, nearest_position) # tuple of nearest array position and distance to protein	






# list of crispr arrs for each genome
for arr in crispr_arrs:
	arr_entry = arr.id.split("[")
	arr_entry = arr_entry[0]
	if (" " in arr_entry):
		arr_entry.split(" ")
		arr_entry = arr_entry[0] # if no space then entire string
=	if (arr_entry in crispr_arr_dict):
		crispr_arr_dict[arr_entry].append(arr)
	else:
		crispr_arr_dict[arr_entry] = [arr]
index_pattern = re.compile("::[0-9]*:[0-9]*")

ret_seq = []

i = 0
for sequ in sequences:
	returned_seq = re.search(index_pattern, sequ.description)

	genome_id = returned_seq.group(0)
	contig_start_coordinate = genome_id
=	genome_id = sequ.description.split(genome_id)

	genome_id = genome_id[0]
	genome_id = genome_id.split("|")
	genome_id = genome_id[-1]
	contig_start_coordinate = contig_start_coordinate[2:]
	contig_start_coordinate = contig_start_coordinate.split(":")
	contig_start_coordinate = contig_start_coordinate[0] # first coordinate of contig
	
	if (genome_id == None):
		print("error in pattern matching!")
	genome_id = str(genome_id)
	if (genome_id in crispr_arr_dict):
		crispr_arr_sequences = crispr_arr_dict[genome_id] # id for crispr array including position
		# make distance calculation here! 
		protein_position = sequ.id.split("|")
		protein_position2 = protein_position [5]
		protein_position = protein_position[4]
		protein_position = int(protein_position)
		protein_position = (protein_position + int(protein_position2)) / 2

		arr_positions = []
		for arr_entry in crispr_arr_sequences:


			arr_entry = arr_entry.description.split("Pos=")
			arr_entry = arr_entry[1]
			arr_entry = arr_entry.split("]")
			arr_entry = arr_entry[0]
			arr_positions.append(int(arr_entry))

		nearest_arr = closest_crispr_array(arr_positions, protein_position, int(contig_start_coordinate))
		nearest_arr_position = nearest_arr[0]
		nearest_arr_distance = nearest_arr[1]


		sequ_cpy = sequ
		sequ_cpy.description += " " + "CRISPR_Position=" + str(nearest_arr_position) + " " + "Distance=" + str(nearest_arr_distance) 
		ret_seq.append(sequ_cpy)
		i += 1
		if (i % 1000 == 0):
			print(i)
	else:	
		print("error, can't locate array from id!")
		print(genome_id, sequ.description)

SeqIO.write(ret_seq, "position_n_distance_" + sys.argv[1], "fasta")