# filtered_arrs_appender

# script to do 3 things:
# 1. Number by CRISPR-array. Assign spacer positions relative to the PFS.
# 2. Retrieve only the rows with corresponding filtered hits (but keep the starting mapped hits)
# 3. Index all possible genome_id_spacer_entries to the Genome_id and spacer. -> simplest is to map by genome_id-spacer coords. If there is a conflict then use the nearest spacer coordinates. There shouldn't ben one though!!! BEcause the start positions are drawn directly from the spacer table.

import sys
import csv
import copy

# iteratively determine if the start and end coordinates overlap, then make these the new bounded coordinates
def duplicates(arr_mode_seq,stend):
	first_seqs = {}
	duplicates = {}
	index = 1
	if (stend == "start"):
		index = 0
	for ele in arr_mode_seq:
		if ele[index] not in first_seqs:
			first_seqs[ele[index]] = ele
		else:
			if ele[index] not in duplicates:
				duplicates[ele[index]] = ele
	return duplicates	

def correct(contig_id): # This may be the dictionary directly
	# iterate over the set rows.
	# If either the start or end coordinate is the same then assign the array start and end coordinates as the lowest and highest values possible from the set of values in the dataset
	# This is important to ensure that the final (PFS) spacer is not munted.
	hit_table = contig_id
	contig_dict = {}
	for row in hit_table:
		if row[0] not in contig_dict:
			contig_dict[row[0]] = [row]
		else:
			contig_dict[row[0]].append(row)
	new_arrs = []
	for array in contig_dict:
		
		forward = 0 # Left = forward, Right = Reverse 
		reverse = 0
		start_mode = []
		arr_end = []
		for row in copy.deepcopy(contig_dict[array]):
			start_mode.append((row[6],row[7]))
			
		arr_start_mode = duplicates(start_mode, "start")
		arr_end_mode = duplicates(start_mode,"end")
		for row in contig_dict[array]:
			if ((row[6] not in arr_start_mode and row[7] in arr_end_mode)):
		#		if (row[7] == "271"):

				row[6] = str(arr_end_mode[row[7]][0])
			new_arrs.append(row)

	return new_arrs

def overlap (row, dict_keys):
	start = int(row[7])
	end = int(row[8])
	
	for coord in dict_keys:
		coord_start = int(coord[0])
		coord_end = int(coord[1])
		if (coord_start < start < coord_end or coord_start < end < coord_end):
			if (start < coord_start):
				new_coord_start = start 
			else:
				new_coord_start = coord_start
			if (end < coord_end):
				new_coord_end = end 
			else:
				new_coord_end = coord_end  
			return (start, end, new_coord_start, new_coord_end, coord_start, coord_end, row)	
		
	return (start, end, "NA", "NA", coord_start, coord_end, row)
							

def is_overlapping (crispr_array):
	# only one row found!!
	start_outer_bound = int(crispr_array[0][7]) # check this!!
	end_outer_bound = int(crispr_array[0][8])
	ret_dict = {}
	ret_dict[(int(start_outer_bound), int(end_outer_bound))] = crispr_array[0][0] # check this!!
	for row in crispr_array:
		start = int(row[6])
		end = int(row[7])
		new_key_value_pair = overlap(row, ret_dict.keys())

		if (new_key_value_pair[2] != "NA"):
			if (start == new_key_value_pair[2] and end == new_key_value_pair[3]):
				ret_dict[(int(start), int(end))].append(row)
			else:
				existing_value = ret_dict.pop((new_key_value_pair[4], new_key_value_pair[5]))
				ret_dict[(int(new_key_value_pair[2]), int(new_key_value_pair[3]))] = existing_value
				ret_dict[(int(new_key_value_pair[2]), int(new_key_value_pair[3]))].append(row)
		else:
			ret_dict[(int(start), int(end))] =  [row] 
	return list(ret_dict.values())		
			


# would be safer to use ordered_dict (as dict is only ordered in python 3.7 or later)
def extract_spacer_coords(array):
	ret_list = []
	for row in array:
		ret_list.append((row[10],row[11]))
	return ret_list

# generate the coordinates of the largest concensous array
def anneal(arr_keys):
	contig_coords = []
	for arr_coord in arr_keys:
		index = 0
		# for loop can't store data!!!
		a = 0
		while (a < len( contig_coords)):
			contig_coords_start = int(float(contig_coords[a][0][0]))
			contig_coords_end = int(float(contig_coords[a][0][1]))

		#	if (int(float(arr_coord[0])) <= contig_coords_start <= int(float(arr_coord[1])) or int(float(arr_coord[0])) <= contig_coords_end <= int(float(arr_coord[1])) or (int(float(arr_coord[0])) <= contig_coords_start and int(float(arr_coord[1]) >= contig_coords_end))):
			if (contig_coords_start <= int(float(arr_coord[0])) <= contig_coords_end or contig_coords_start <= int(float(arr_coord[1])) <= contig_coords_end or (contig_coords_start <= int(float(arr_coord[0])) and contig_coords_end >= (int(float(arr_coord[1])))) or (contig_coords_start >= int(float(arr_coord[0])) and contig_coords_end <= (int(float(arr_coord[1]))))):

				if (int(float(arr_coord[0])) <= contig_coords_start):
					
					contig_coords[a][0][0] = int(float(arr_coord[0]))
				if (int(float(arr_coord[1])) >= contig_coords_end):	
					contig_coords[a][0][1] = int(float(arr_coord[1]))

				contig_coords[a][1].append(arr_coord)
				index = 1
			a += 1		
		if (index == 0):	
			contig_coords.append([list(arr_coord),[arr_coord]])

	return contig_coords

# create a single array as the largest non-overlapping collection of spacers from CRISPR-CRT, PILER-CR and CRISPRdetect
def largest_array_merge(hit_table):
	contig_dict = {}

	for row in hit_table:
		if row[0] not in contig_dict:
			contig_dict[row[0]] = [row]
		else:
			contig_dict[row[0]].append(row)
	new_arrs = []
	# run through each contig and find the set of array start and end coordinates
	for array in contig_dict:
		# create an array dict nested inside the contig dict
		# If the start/end array coordinates overlap, then merge the two array dictionaries
		array_dict = {} # dict for each individual array in the contig
		
		for row in contig_dict[array]:
			if (row[6], row[7]) not in array_dict:
				array_dict[(row[6], row[7])] = [row]
			else:
				array_dict[(row[6], row[7])].append(row)
		arr_keys = list(array_dict.keys())
	#	overlapping_array_pairs = overlaps(arr_keys) # This is a valid if inefficent way to access the array keys.
		# add new array coordinates
		annealled_groups = anneal(arr_keys)
		for group in annealled_groups:
			start_coord = str(group[0][0])
			end_coord = str(group[0][1])
			for problem_arr in group[1]:
				i = 0
				while (i < len(array_dict[problem_arr])):
					array_dict[problem_arr][i][6] = start_coord
					array_dict[problem_arr][i][7] = end_coord
					i += 1

		for arr in array_dict:
			for row in array_dict[arr]:
				new_arrs.append(row)
	
	return new_arrs

# Get rid of any spacers duplicated by prediction with multiple array prediction tools by checking for overlaps.
def is_in (row, array, offset, arr_dict):
	row_start = int(float(row[10]))
	row_end = int(float(row[11]))
	for arr_row in array:
		# not sure whether == is powerful enough in python to perform this operation on it's own. May need custom function
		if (row == arr_row):
			print("skip!!")
			continue
		# must not include the self mapping coord
		if (row_start <= int(float(arr_row[10])) <= row_start + offset or row_start - offset <= int(float(arr_row[10])) <= row_start or  row_end <= int(float(arr_row[11])) <= row_end + offset or row_end - offset <= int(float(arr_row[11])) <= row_end):
			return True
	return False		



def lam_sort(x):
#	print(x)
	return int(float(x[10]))

# Number spacers from the PPS
def number_arrays (merged_table): # number arrays can assign the relative position of each spacer
	ret_table = []
	array_dict = {}
	arr_index = 0
# 1. First create a dictionary with the array start/end as the keys and spacer rows as the values.
	for row in merged_table:
		genome_id = row[0].split(" ") [0]
		if (genome_id + row[6] + row[7] not in array_dict):
			array_dict[genome_id + row[6] + row[7]] = [row]# want to be able to add different arrays in each genome as seperate groups
			arr_index += 1
		else: # This is the problematic row as with arr_index incrementing there can never be two identical rows
			# Is this condition ever triggered?
		#	print("Yello!!")
			array_dict[genome_id + row[6] + row[7]].append(row)
		#	print(array_dict[genome_id + row[6] +row[7]])
	i = 0 
	# iterate through the array containing a list of lists of spacer coordinates for each CRISPR_array. Add a row if the spacer value is not within within 10bp of the other rows. Must not map to self  
# 2. NOT WHAT THIS DOES!!! - (Determine if the spacer start/end coordinates are all within the bounds of one array!!)
# consider whether this function is worth it? May just unnessessarily distort data	
		
	#	Is this function nessessary. A simplier way to create an array dictionary would be to just use the start/end coordinates. Have already previously checked for overlap!!
	#	array_dict[array] = is_overlapping(array_dict[array]) # should split list entries in seperate arrays in individual genomes. Need to test this!!!
		#

    # 3. Sort the arrays and number by the oldest spacers!!
	crispr_index = 0
	for genome in array_dict:
		if (array_dict[genome][0][1] == "Forward"):
			try:
				array_dict[genome] = sorted(array_dict[genome], key=lam_sort, reverse=True) # sort by descending spacers because PFS is the last spacer
			except:
				print(array_dict[genome])
		# Otherwise must be reverse
		else:
		#	print(array_dict[genome])
			array_dict[genome] = sorted(array_dict[genome], key=lam_sort)	
		k = 0 # This IS THE SPACER NUMBER INDEX IN THE PRESUMED ORDER OF WHEN SPACERS WERE INTEGRATED BASED ON THEIR ORIENTATION!! IS THIS A VALID ASSUMPTION???
		while(k < len(array_dict[genome])):
			array_dict[genome][k].append(crispr_index) # append the crispr_array number
			array_dict[genome][k].append(str(k))
			k += 1
		crispr_index += 1
		ret_table.append(array_dict[genome])


	return ret_table # may want to save this to seperate file!!!

# function to generate a concensous set of array predictions and to number each spacer from the PPS
def spacer_hits_appender(filtered_spacer_hit_url, merged_spacer_table_url, output_url):
	filtered_spacer_hit_url = open(filtered_spacer_hit_url, "r")
	merged_spacer_table_url = open(merged_spacer_table_url, "r")
	output_url = open(output_url, "w")

	
	merged_spacer_table = list(csv.reader(merged_spacer_table_url))
	header = merged_spacer_table[0]
	# unsure if filtering these out is nessessarily the best approach? 
	i = 0
	while (i < len(merged_spacer_table)):
		if (merged_spacer_table[i][10] == ''):
			merged_spacer_table.pop(i)
		else:
			i += 1	

	merged_spacer_table = largest_array_merge(merged_spacer_table)		
	merged_spacer_table = number_arrays(merged_spacer_table[1:]) # assign numbers for each spacer relative to the PFS
	filtered_spacer_hits = list(csv.reader(filtered_spacer_hit_url))
	'''
	print(merged_spacer_table)
	out_file = open("merge_numbered.csv","w")
	number_writer = csv.writer(out_file)
	for ele in merged_spacer_table:
		for row in ele:
			number_writer.writerow(row)
	out_file.close() 
#	exit()
	ret_out = open("check_merged_numbering.csv","w")
	spoof_writer = csv.writer(ret_out)
	for crispr in merged_spacer_table:
		for row in crispr:
			spoof_writer.writerow(row)
	ret_out.close()
'''

	spam_writer = csv.writer(output_url)
	spam_writer.writerow(["Spacer_id","Phage_id","Perc_id","Length", "Mismatches","Gapopen","query_start","query_end","Mapped_start_site","Mapped_end_site","evalue","bitscore","Genome_id","orientation","orientation_score","orientation_confidence","questionable_array","array_score","CRISPR-start","CRISPR-end","repeat_start","repeat_end","spacer_start","spacer_end","dr_repeat_original","dr_repeat_concensous","spacer","Array_tool","run","array_number","spacer_number"])
	out_table = []
	out_row = []

	spacer_dict = {}
	# why is this double looped!!
	for hit in merged_spacer_table:
		for spacer in hit:
			spacer_id = spacer[0].split(" ")
			spacer_id = spacer_id[0]
			spacer_dict[(spacer_id.strip(), str(int(float(spacer[10]))).strip(), str(int(float(spacer[11]))).strip())] = spacer
#	all_keys = spacer_dict.keys()
#	my_out = open("outfile.csv","w")
#	spa_writer = csv.writer(my_out)
#	for akey in all_keys:
#		spa_writer.writerow(akey)
#	my_out.close()

	i = 0
	for spacer in filtered_spacer_hits:
	#	print("start of loop")
		spacer_id = spacer[0].split("|")
	#	print(spacer_id)
		spacer_start = spacer_id[1].split("spacer_start_pos:") [1]
		spacer_end = spacer_id[2].split("spacer_end_pos:") [1]
		spacer_id = spacer_id[0]

		spacer_start = spacer_start.split("|") [0]
		spacer_end = spacer_end.split("|") [0]
		if (spacer_id.strip(), spacer_start.strip(), spacer_end.strip()) in spacer_dict:
			# This might be a flawed line!!
			out_row = spacer + spacer_dict[(spacer_id.strip(), spacer_start.strip(), spacer_end.strip())]
		#	print("Checked_row!!")
		else:
			# this needs to be subject to a check as well!!
			print("Error, original mapping sequence and orientation not found") # unless a mistake has been made, there should never be a case where the spacer + orientation information can't be retrieved and mapped to the hits. However, some sequences are removed from the merged table. Which may account for this!!
			print(spacer_id.strip(), spacer_start.strip(), spacer_end.strip())
			print("end of id")
		spam_writer.writerow(out_row)
		out_table.append(out_row)
		i += 1
	#	print(len(filtered_spacer_hits))
	#	print(i)
	# First number the arrays (may want to save to seperate file) then merge with filtered spacer hits -> Use existing code as a template!!
	print("done_arr_appending!!")
	filtered_spacer_hit_url.close()
	merged_spacer_table_url.close()
	output_url.close()
	return 0