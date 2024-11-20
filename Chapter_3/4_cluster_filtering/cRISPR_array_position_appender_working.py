# CRISPR array position appender

# the plan.

# load relabelled file
# load concatenated CRISPR-array file
# extract position score
# Find the minimum distance between the array and the starting point of the arrays
# this requires computing the abolsute positions of the proteins from the CRISPR array
# the abolsute position will be the position given by prodigal/genemark + starting array position
# if strand is negative the starting position should have a higer score than the negative position
#multiple mapping between arrays and ids mean that the minimal distance between the ids and arrays must be found

from Bio import SeqIO
import csv
import sys
import re


sequences = SeqIO.parse(sys.argv[1], "fasta")
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
	#print(arr_entry)
	arr_entry = arr_entry[0]
	if (" " in arr_entry):
		arr_entry.split(" ")
	#	print(arr_entry)
		arr_entry = arr_entry[0] # if no space then entire string
#		print(arr_entry)
	if (arr_entry in crispr_arr_dict):
		crispr_arr_dict[arr_entry].append(arr)
		#print("old entry")
	else:
		crispr_arr_dict[arr_entry] = [arr]
		#print("new_entry")
#need to check this pattern
index_pattern = re.compile("::[0-9]*:[0-9]*")

ret_seq = []

# 
i = 0
for sequ in sequences:
	returned_seq = re.search(index_pattern, sequ.description)
#	print(sequ.description)

#	print("regex working?")
	genome_id = returned_seq.group(0)
#	print(genome_id)
#	print(genome_id)
	contig_start_coordinate = genome_id
#	print(contig_start_coordinate)
	genome_id = sequ.description.split(genome_id)

	genome_id = genome_id[0]
	genome_id = genome_id.split("|")
	genome_id = genome_id[-1]
	#print(genome_id, sequ.description)
	contig_start_coordinate = contig_start_coordinate[2:]
	contig_start_coordinate = contig_start_coordinate.split(":")
#	print("contig start coordinate:")
#	print(contig_start_coordinate)
	contig_start_coordinate = contig_start_coordinate[0] # first coordinate of contig
	#genome_id = genome_id[1:] # genome id
	
	if (genome_id == None):
		print("error in pattern matching!")
	genome_id = str(genome_id)
	if (genome_id in crispr_arr_dict):
		crispr_arr_sequences = crispr_arr_dict[genome_id] # id for crispr array including position
	#	print("mapped!")
		# make distance calculation here! 
		# should be 
		# need to find absolute position of protein
		protein_position = sequ.id.split("|")
		protein_position2 = protein_position [5]
		protein_position = protein_position[4]
		protein_position = int(protein_position)
		protein_position = (protein_position + int(protein_position2)) / 2

		arr_positions = []
		for arr_entry in crispr_arr_sequences:


			arr_entry = arr_entry.description.split("Pos=")
		#	print (crispr_arr_sequences)
		#	print(arr_entry)
			arr_entry = arr_entry[1]
		#	print(arr_entry)
			arr_entry = arr_entry.split("]")
			arr_entry = arr_entry[0]
			arr_positions.append(int(arr_entry))

		nearest_arr = closest_crispr_array(arr_positions, protein_position, int(contig_start_coordinate))
		nearest_arr_position = nearest_arr[0]
		nearest_arr_distance = nearest_arr[1]


		sequ_cpy = sequ
		sequ_cpy.description += " " + "CRISPR_Position=" + str(nearest_arr_position) + " " + "Distance=" + str(nearest_arr_distance) 
	#	print (sequ_cpy)
		ret_seq.append(sequ_cpy)
		i += 1
		if (i % 1000 == 0):
			print(i)
		#print("Good!")
	else:	
		print("error, can't locate array from id!")
		print(genome_id, sequ.description)

SeqIO.write(ret_seq, "position_n_distance_" + sys.argv[1], "fasta")