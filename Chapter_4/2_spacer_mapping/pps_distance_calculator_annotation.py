# spacer PPS distance calculator


import sys
import csv
import copy

def forward_sorter(row):
	return (float(int (row[20]) + int (row[21]) / 2))

# identify the PPS, then compute the distance from the PPS to the spacers. Spacers should be assigned based on whether spacers map to either the forward or reverse strand.
# note, an updated version of this code was used for spacer distribution analysis. This script is included for backward-compatibility with the annotation workflow.
def pps_compute(input_url):
	print("Start:")
	with open (input_url, "r") as csvfile:
		final_mapping_table = list(csv.reader(csvfile))

	# first filter NA rows 
	ret_out = open(input_url + "_no_NA.csv", "w")
	spam_writer = csv.writer(ret_out)
	for row in final_mapping_table:
		if (row[8] != "NA" and row[8] != ''): # may need to change this row number
			spam_writer.writerow(row)
	ret_out.close()



	with open (input_url + "_no_NA.csv", "r") as csvfile:
		final_mapping_table = list(csv.reader(csvfile))

	spacer_host_hits_dict = {}
	header = final_mapping_table[0]
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

	spacer_host_hits = list(spacer_host_hits_dict.values())
	# now spacer-host pairs are grouped.
	# Need to go through each row. Select the PPS. Compute the distance from the PPS for each entry and assign whether or not the hit is on the same, or complementary strand.
	ret_out = open(input_url + "_distances_annotated.csv", "w")
	spam_writer = csv.writer(ret_out)
	spam_writer.writerow(["Spacer_id","Phage_id","Perc_id","Length", "Mismatches","Gapopen","query_start","query_end","Mapped_start_site","Mapped_end_site","evalue","bitscore","Genome_id","orientation","orientation_score","orientation_confidence","questionable_array","array_score","CRISPR-start","CRISPR-end","repeat_start","repeat_end","spacer_start","spacer_end","dr_repeat_original","dr_repeat_concensous","spacer","Array_tool","RUN","array_number","spacer_number","distance","mapped_strand"])
	print("Hi!!")
	# add header!!
	for spacer_host_hit in spacer_host_hits:
		array_orientation = spacer_host_hit[0][13]
		if (array_orientation == "Forward"):
			pps = list(sorted(spacer_host_hit, reverse=True,key=forward_sorter)) [0] # make sure that the sequences are not already sorted!!
			# THESE INDEXES WILL NEED TO BE CHANGED - Actually this works out!!
			spacer_start = pps[8]
			spacer_end = pps[9]
			pps_coord = float((int(spacer_start) + int(spacer_end)) / 2)
			for hit in spacer_host_hit:
				# SPACER_TABLE_INDEX WILL NEED TO BE MODIFIED
				spacer_start_hit = hit[8]
				spacer_end_hit = hit[9]
				midpoint_spacer_distance = float ((int(spacer_start_hit) + int(spacer_end_hit)) / 2)
				midpoint_spacer_distance = midpoint_spacer_distance - pps_coord
				
				if (int(spacer_start_hit) <= int(spacer_end_hit)):
					strand = "1"
				else:
					strand = "-1"
				my_hit = copy.deepcopy(hit)
				my_hit.extend([str(midpoint_spacer_distance), strand])
				spam_writer.writerow(my_hit)
		else: # spacer_host_hit == "Reverse"
		#	print(spacer_host_hit)
			pps = list(sorted(spacer_host_hit, reverse=False,key=forward_sorter)) [0]
			pps_default = list(sorted(spacer_host_hit, reverse=True,key=forward_sorter)) [0]
			# THESE INDEXES WILL NEED TO BE CHANGED 

			spacer_start = pps[8]
			spacer_end = pps[9]

			pps_coord = float ((int(spacer_start) + int(spacer_end)) / 2)
			for hit in spacer_host_hit:
				# SPACER_TABLE_INDEX WILL NEED TO BE MODIFIED!!
				spacer_start_hit = hit[8]
				spacer_end_hit = hit[9]

				midpoint_spacer_distance = float ((int(spacer_start_hit) + int(spacer_end_hit)) / 2)
				midpoint_spacer_distance = midpoint_spacer_distance - pps_coord
				if (int(spacer_start_hit) <= int(spacer_end_hit)):
					strand = "1"
				else:
					strand = "-1"
				my_hit = copy.deepcopy(hit)
				my_hit.extend([str(midpoint_spacer_distance), strand])
				spam_writer.writerow(my_hit)
	ret_out.close()
	return 0
