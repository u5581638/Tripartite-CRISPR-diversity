# spacer PPS distance calculator
# identify the PPS, then compute the distance from the PPS to the spacers.


import sys
import csv
import copy

def forward_sorter(row):
	return (float(float (row[20]) + float (row[21]) / 2))


# identify the PPS, then compute the distance from the PPS to the spacers. Spacers are assigned based on whether spacers map to either the forward or reverse strand.
def pps_compute(input_url):
	print("Start:")
	csvfile = open (input_url, "r") 
	final_mapping_table = csv.reader(csvfile)
	# first filter NA rows 
	ret_out = open(input_url + "_no_NA.csv", "w")
	spam_writer = csv.writer(ret_out)
	spam_writer.writerow(["Spacer_id","Phage_id","Perc_id","Length", "Mismatches","Gapopen","query_start","query_end","Mapped_start_site","Mapped_end_site","evalue","bitscore","Genome_id","orientation","orientation_score","orientation_confidence","questionable_array","array_score","CRISPR-start","CRISPR-end","repeat_start","repeat_end","spacer_start","spacer_end","dr_repeat_original","dr_repeat_concensous","spacer","Array_tool","RUN","array_number","spacer_number","distance","mapped_strand"])
	for row in final_mapping_table:
		if (row[8] != "NA" and row[8] != ''): # may need to change this row number
			spam_writer.writerow(row)
	csvfile.close()
	ret_out.close()

	print("NO_NA written to drive")

	my_in = open (input_url + "_no_NA.csv", "r")
	final_mapping_table = csv.reader(my_in)
	print("Dict start:")
	spacer_host_hits_dict = {}
	header = next(final_mapping_table)
	final_mapping_table = final_mapping_table
	for row in final_mapping_table:
		# THESE INDEXES WILL NEED TO BE REPLACED WITH THE UPDATED INDEXES IN THE REAL TABLE!!
		genome_id = row[0].split("|") [0]
		mapped_hit = row[1]
		crispr_start = row[18] 
		crispr_end = row[19]
		# END OF MODIFICATIONS!!

		spacer_host_id = genome_id.strip() + " " + mapped_hit.strip() + " " + crispr_start.strip() + " " + crispr_end.strip()
		if (spacer_host_id not in spacer_host_hits_dict):
			spacer_host_hits_dict[spacer_host_id] = [row]
		else:
			spacer_host_hits_dict[spacer_host_id].append(row)
	my_in.close()
	spacer_host_hits = list(spacer_host_hits_dict.values())
	# now spacer-host pairs are grouped.
	# Need to go through each row. Select the PPS. Compute the distance from the PPS for each entry and assign whether or not the hit is on the same, or complementary strand.
	print("Final_IO")
	ret_out = open(input_url + "_distances_annotated.csv", "w")
	spam_writer = csv.writer(ret_out)
	spam_writer.writerow(["Spacer_id","Phage_id","Perc_id","Length", "Mismatches","Gapopen","query_start","query_end","Mapped_start_site","Mapped_end_site","evalue","bitscore","Genome_id","orientation","orientation_score","orientation_confidence","questionable_array","array_score","CRISPR-start","CRISPR-end","repeat_start","repeat_end","spacer_start","spacer_end","dr_repeat_original","dr_repeat_concensous","spacer","Array_tool","RUN","array_number","spacer_number","distance","mapped_strand"])
	# add header!!
	for spacer_host_hit in spacer_host_hits:
		array_orientation = spacer_host_hit[0][13]
		if (array_orientation == "Forward"):
			pps = list(sorted(spacer_host_hit, reverse=True,key=forward_sorter)) [0] # make sure that the sequences are not already sorted!!
			spacer_start = pps[8]
			spacer_end = pps[9]
			pps_coord = float((int(spacer_start) + int(spacer_end)) / 2)
			pps_sense = 1
			if (int(spacer_start) > int(spacer_end)):
				pps_sense = -1
			for hit in spacer_host_hit:
				# SPACER_TABLE_INDEX WILL NEED TO BE MODIFIED
				spacer_start_hit = hit[8]
				spacer_end_hit = hit[9]
				target_sense = 1 
				if (int(spacer_end_hit) < int(spacer_start_hit)):
					target_sense = -1
				midpoint_spacer_distance = float ((int(spacer_start_hit) + int(spacer_end_hit)) / 2) # could be as simple as adding abs here?
				if (pps_sense == 1 and target_sense == 1):
					midpoint_spacer_distance = midpoint_spacer_distance - pps_coord
					strand = "1"
				elif (pps_sense == 1 and target_sense == -1):
					midpoint_spacer_distance = midpoint_spacer_distance - pps_coord
					strand = "-1"	
				elif (pps_sense == -1 and target_sense == 1):
					midpoint_spacer_distance = pps_coord - midpoint_spacer_distance
					strand = "-1"
				elif (pps_sense == -1 and target_sense == -1):
					midpoint_spacer_distance = pps_coord - midpoint_spacer_distance			
					strand = "1"
				else:
					print("Error!!")

				my_hit = copy.deepcopy(hit)
				my_hit.extend([str(midpoint_spacer_distance), strand])
				spam_writer.writerow(my_hit)
		else: # spacer_host_hit == "Reverse"
			pps = list(sorted(spacer_host_hit, reverse=False,key=forward_sorter)) [0]
			pps_default = list(sorted(spacer_host_hit, reverse=True,key=forward_sorter)) [0]
			spacer_start = pps[8]
			spacer_end = pps[9]
			pps_coord = float ((int(spacer_start) + int(spacer_end)) / 2)
			pps_sense = -1

			if (int(spacer_start) > int(spacer_end)):
				pps_sense = 1
			for hit in spacer_host_hit:
				# SPACER_TABLE_INDEX WILL NEED TO BE MODIFIED!!
				spacer_start_hit = hit[8]
				spacer_end_hit = hit[9]
				target_sense = -1 
				if (int(spacer_end_hit) < int(spacer_start_hit)):
					target_sense = 1
				midpoint_spacer_distance = float ((int(spacer_start_hit) + int(spacer_end_hit)) / 2)
				# compute midpoint distance
				if (pps_sense == 1 and target_sense == 1):
					midpoint_spacer_distance =  midpoint_spacer_distance - pps_coord
					strand = "1"
				elif (pps_sense == 1 and target_sense == -1):
					midpoint_spacer_distance = midpoint_spacer_distance - pps_coord
					strand = "-1"
				elif (pps_sense == -1 and target_sense == 1):
					midpoint_spacer_distance = pps_coord - midpoint_spacer_distance
					strand = "-1"

				elif (pps_sense == -1 and target_sense == -1):
					midpoint_spacer_distance = pps_coord - midpoint_spacer_distance		
					strand = "1"
				else:
					print("Error!!")
				my_hit = copy.deepcopy(hit)
				my_hit.extend([str(midpoint_spacer_distance), strand])
				spam_writer.writerow(my_hit)
	ret_out.close()
	return 0
pps_compute(sys.argv[1])

# INPUT: hitmap table containing only spacer matches to 2 or more sites in target contigs (from the same array)
# i.e. cas12a.fasta_all_hits.csv_genomes.fasta_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_hitmap.csv_standardised.csv_non_redundant.csv_filtered.csv_2_or_more_hits.csv
# OUTPUT: hitmap table annotated each hits appending the distance from the PPS. The PPS itself is removed from the table. Strand directionality with respect to the PPS is also given.
# i.e. cas12a.fasta_all_hits.csv_genomes.fasta_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_hitmap.csv_standardised.csv_non_redundant.csv_filtered.csv_2_or_more_hits.csv_distances_annotated.csv
# SHELL: python3 pps_distance_calculator_annotation_corrected_strand_direction.py cas12a.fasta_all_hits.csv_genomes.fasta_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_hitmap.csv_standardised.csv_non_redundant.csv_filtered.csv_2_or_more_hits.csv