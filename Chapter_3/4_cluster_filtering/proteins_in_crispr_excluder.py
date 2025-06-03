# proteins_in_crispr_excluder.py

from Bio import SeqIO 
import sys
import csv

# plan:

# first need to load the sequences and positions.
# positions should be in dictionary form. May need to proccess the keys to be recognisable via sequences
# lookup via sequence should be sufficent
# then need to check for each sequence that the crispr arr does not span any of the coordinates given by pilercr

# INPUT: list of protein representatives from each cluster with start and end sites labelled in header. Also requires a table containing the contig ides and their CRISPR-array coordinates.
# i.e. 1. representative_distance_5000_10000avg_dist_position_n_distance_relabelled_grt_3_grt_3_cluster_declustered_one_member_gt_300_renumbered_sorted_descending_deduplicated_flattened_all_sequences.fasta_p_formatted.fasta 
# 	   2. spliced_cleaned_rerun_20kb_local_windows_combined.fasta_pilercr.fa_real_arr_positions.csv
# OUTPUT: Set of protein representatives which are not CRISPR-array ORFs
# i.e. representative_distance_5000_10000avg_dist_position_n_distance_relabelled_grt_3_grt_3_cluster_declustered_one_member_gt_300_renumbered_sorted_descending_deduplicated_flattened_all_sequences.fasta_p_formatted.fastanot_in_crisprs.fasta
# SHELL: python3 proteins_in_crispr_excluder.py representative_distance_5000_10000avg_dist_position_n_distance_relabelled_grt_3_grt_3_cluster_declustered_one_member_gt_300_renumbered_sorted_descending_deduplicated_flattened_all_sequences.fasta_p_formatted.fasta spliced_cleaned_rerun_20kb_local_windows_combined.fasta_pilercr.fa_real_arr_positions.csv
# import the sequences

sequences = SeqIO.parse(sys.argv[1], "fasta")
with open(sys.argv[2], "r") as csvfile:
	crispr_coords = list(csv.reader(csvfile))
ret_list = []
crispr_coords_dict = {}
# dicionary of crispr records
print(len(crispr_coords))
for record in crispr_coords:
	record_id = record[0].split("|")
	record_id = record_id[-1]
	if record_id not in crispr_coords:
		crispr_coords_dict[record_id] = [record]
	else:
		crispr_coords_dict[record_id].append(record)

print(crispr_coords_dict.keys())
for sequence in sequences:
# extract the protein start and end sites
	sequence_description = sequence.id.split("|")
	sequence_id = sequence_description[0]
	sequence_start_pos = int(sequence_description[4])
	sequence_end_pos = int(sequence_description[5])
	if (sequence_id in crispr_coords_dict):
		arr_coords = crispr_coords_dict[sequence_id]
		for arr_coord in arr_coords: # arr_coord[0]: arr start pos. arr_coord[1]: arr end pos
			if (sequence_start_pos < int(arr_coord[1]) and sequence_end_pos < int(arr_coord[1])):
				ret_list.append(sequence)
			elif (sequence_start_pos > int(arr_coord[2]) and sequence_end_pos > int(arr_coord [2])):
				ret_list.append(sequence)
			else: # arr sequence somewhere within protein
				break		

	else:
		print(sequence_id + " is missing from arrs")	
SeqIO.write(ret_list, sys.argv[1] + "not_in_crisprs.fasta", "fasta")		