# program to reconcile CRISPR-CRT, PILERCR and CRISPRdetect predictions into a single table

import csv
import sys
import copy

# Test this program by running function by function.
# Use a small sample of 10 individual genomes/contigs

# Note that dicts should be ordered dicts


def add_array (crt_array):
	ret_list = []
	dr_con = []
	for array in crt_array:

		ret_list.append([array[0],array[13], "NA", "NA","NA","NA", array[2], array[3], array[4], array[5], array[6], array[7], array[8], "Not_determined", array[9], "CRT"])
		dr_con.append(array[8])
	con_dr = max(set(dr_con), key=dr_con.count)
	i = 0
	while (i < len(ret_list)):
		ret_list[i][13] = con_dr
		i += 1
	return ret_list	


def nearest_array (crt_spacer_start, crt_spacer_end, array):
	distances = {}
	for row in array:
		crispr_start = int(row[6])
		crispr_end = int(row[7])
		distance = abs(abs((crt_spacer_start + crt_spacer_end)/2) - abs((crispr_start + crispr_end)/2))
		distances[distance] = row
	shortest_distance = min(distances.keys())
	return distances[shortest_distance]	

# updated implementation of construct_crisprdetect_row_minimal
def construct_crisprdetect_row_minimal_update (row, nearest_row,crisprdetect_dict_array):
	ret_row = []
	repeat_start = int(row[4])
	repeat_end = int(row[5])
#	start_distance = int(repeat_start) - int(nearest_row[6])
# 	end_distance = int(repeat_end) - int(nearest_row[7]) 
	# extend start and end array coordinates 
	if (nearest_row[6] < row[2] < nearest_row[7] and nearest_row[6] < row[3] < nearest_row[7]):
		crispr_start = nearest_row[6]
		crispr_end = nearest_row[7]
	#	print("Mate!!")
	elif (nearest_row[6] < row[2] < nearest_row[7]):
		crispr_start = nearest_row[6]
		crispr_end = row[3]
	#	print("Interesting")
		i = 0
		while (i < len(crisprdetect_dict_array)): # only update entries wher crisprdetect_dict_array = nearest row
			if (crisprdetect_dict_array[i][7] == nearest_row[7] ):
				crisprdetect_dict_array[i][7] = crispr_end
			i += 1
	elif (nearest_row[6] < row[3] < nearest_row[7]):
	#	print("Sure!!")
		crispr_start = row[2]
		crispr_end = nearest_row[7]
		i = 0
		while (i < len(crisprdetect_dict_array)): # only update entries wher crisprdetect_dict_array = nearest row
			if (crisprdetect_dict_array[i][6] == nearest_row[6] ):
				crisprdetect_dict_array[i][6] = crispr_start
			i += 1
	else:
		# if either the crispr_start or crispr_end coordinates match in crispr_detect_dict then update these rows
		crispr_start = row[2]
		crispr_end = row[3]

	ret_row = [[row[0]] + nearest_row[1:6] + [crispr_start] + [crispr_end] + row[4:9] + [nearest_row[13]] + [row[9]] + ["CRT"]]	
	return (ret_row, crisprdetect_dict_array)

# construct a new row with maximal direct repeat corrdinates to be consistent with predictions from CRISPR-CRT
	ret_row = []
	repeat_start = int(row[4])
	repeat_end = int(row[5])
	start_distance = int(repeat_start) - int(nearest_row[6])
	end_distance = int(repeat_end) - int(nearest_row[7]) 
	# extend start and end array coordinates 

	if (row[2] < nearest_row[6]):
		crispr_start = row[2]
	else:
		crispr_start = nearest_row[6]	
	if (row[3] > nearest_row[7]):
		crispr_end = row[3]
	else:
		crispr_end = nearest_row[7]
	ret_row = [[row[0]] + nearest_row[1:6] + [crispr_start] + [crispr_end] + row[4:9] + [nearest_row[13]] + [row[9]] + ["CRT"]]	
	return ret_row




	# compare the start and end array coords. Expand the CRISPR-array range the spacer to be added exceeds the existing bounds.

# problem with this approach is that CRISPR-array self-filtering requires the DR sequence in it's original (not a reconstructed or inferred) form!!
def construct_crisprdetect_row (row, nearest_row):
	ret_row = []
	spacer_start = int(row[6])
	spacer_end = int(row[7])

	start_distance = int(spacer_start) - int(nearest_row[6])
	end_distance = int(spacer_end) - int(nearest_row[7]) 
	stream = min([start_distance, end_distance])

	spacer = row[11]
	repeat = row[10]

	spacer_length = len(spacer)
	repeat_length = len(repeat)
	if (stream == start_distance): # upstream of nearest row
		ret_row = [row[0]] + nearest_row[1:5] + [int(nearest_row[6]) - (spacer_length + repeat_length), nearest_row[7],int(nearest_row[6]) - (spacer_length + repeat_length), int(nearest_row[6]) - (spacer_length), int(nearest_row[6]) - (spacer_length) + 1, int(nearest_row[6]) - 1, row[8], row[9], "CRT"]

	else: # spacer closer to downstream end
		ret_row = [row[0]] + nearest_row[1:5] + [int(nearest_row[6]) , int(nearest_row[7]) + (repeat_length + spacer_length) + 1,int(nearest_row[7]) + 1, int(nearest_row[7]) + repeat_length, int(nearest_row[7]) + repeat_length + 1, int(nearest_row[7]) + 1 + repeat_length + spacer_length, int(nearest_row[6]) - (spacer_length) + 1, row[8], row[11], row[9], "CRT"]
	return ret_row	

# function to reconcile CRISPRdetect and CRISPR-CRT predictions
def crt_reconcile(crisprdetect_table_url, crisprcrt_table_url, output_dir):

	with open(crisprdetect_table_url) as csvfile:
		crisprdetect_table = list(csv.reader(csvfile))
	with open(crisprcrt_table_url) as csvfile2:
		crisprcrt_table = list(csv.reader(csvfile2))

	# 1. Dictionalise all entries into tabular form!!
	# CRISPR_detect table should be pre-inverted!!
	previous_blank_switch = 0
	crisprdetect_dict = {}
	crisprdetect_dict_spacers = {}

	for row in crisprdetect_table: # [1:] should not be nessessary because the inverted table is missing a header
		my_row = copy.deepcopy(row)
	#	print(my_row)
		genome_id = my_row[0].split(" ") [0] # genome_id only

		if genome_id not in crisprdetect_dict:
			crisprdetect_dict[genome_id] = [my_row + ["CRISPRDETECT"]]
			crisprdetect_dict_spacers[genome_id] = {my_row[14]}
		else:
			crisprdetect_dict[genome_id].append(my_row + ["CRISPRDETECT"])
			if (my_row[14] != ''):
				crisprdetect_dict_spacers[genome_id].add(my_row[14])

	crisprcrt_dict = {}

	for row in crisprcrt_table[1:]:
		genome_id = row[0]
		if (genome_id not in crisprcrt_dict):
			crisprcrt_dict[genome_id] = [row] # should this be a dictionary indexes spacers to theur respective rows? No need to create a seperate dictionary doing this mapping and assign spacers as the values of 

		else:
			crisprcrt_dict[genome_id].append(row)

	# Go through this logic with Gaetan/ Tony to check whether it works.
	# crisprcrt_table = copy.deepcopy(crisprcrt_table[1:])
	for array in crisprcrt_dict:
		if (array in crisprdetect_dict):
			# this is actual asking to iterate through the rows of an entry of crisprcrt_dict
			for row in crisprcrt_dict[array]: 
				# Assign CRISPRleader orientation predictions if classification is unconfirmed!!
				if (crisprdetect_dict[array][0][1] == "Unconfirmed"):
					k = 0
					# might need to put in a check here!!
					while (k < len(crisprdetect_dict[array])): 	
						crisprdetect_dict[array][k][1] = row[13] # This shouldn't occur for most spacers is the designation is either forward or reverse
						k += 1

				if ((row[9] in crisprdetect_dict_spacers[array] and row[9] != '') or (previous_blank_switch == 1 and row[9] == '')): # or rev_comp(row[9]) in crisprdetect_dict_spacers[array]): # This should ask if the spacer is in crispr_detect_spacer_set
					# need to append """CRISPRDETECT.CRT""" to row.
					previous_blank_switch = 1
					i = 0
					while i < len (crisprdetect_dict[array]):
						# this needs to be dictionalised to prevent additional iteration!!!!!!!!!!!!!!!!!
						# need to be able to see all rows as list/set in order to be able to decide where to add new spacers. 
						# need a special condition if insertion results in an overlap with another CRISPR-array!!!!
						if (row[9] == crisprdetect_dict[array][i][14]):
							crisprdetect_dict[array][i][15] = "CRISPRDETECT.CRT"
							
						i += 1	

				 # entries are the same!!
				else:
					previous_blank_switch = 0
					# need to iterate through the array to check coordinates!!
					

					crisprcrt_spacer_start = int(row[4])
					crisprcrt_spacer_end = int(row[5])
					i = 0
					while i < len (crisprdetect_dict[array]):
						# this needs to be dictionalised to prevent additional iteration!!!!!!!!!!!!!!!!!
						# need to be able to see all rows as list/set in order to be able to decide where to add new spacers. 
						# need a special condition if insertion results in an overlap with another CRISPR-array!!!!
						crisprdetect_start = int(crisprdetect_dict[array][i][8])
						crisprdetect_end = int(crisprdetect_dict[array][i][9])
						genome_id = crisprdetect_dict[array][i][0].split(" ") [0]
						
						if (crisprdetect_start <= crisprcrt_spacer_start <= crisprdetect_end or crisprdetect_start <= crisprcrt_spacer_end <= crisprdetect_end): # if the spacer is within the array existing array
							# assume this is a pre-existing spacer.
							crisprdetect_dict[array][i][15] = "CRISPRDETECT.CRT" # Add a note that spacers match to both
						else:
					#		print("Hi")

							# need to add to nearest array
							# need to add a CRISPR-array specific number to tables to identify frank arrays of CRT-CRISPRdetect and PILERCR predictions!!!!!!!!!!!!!!
							nearest_row = nearest_array(crisprcrt_spacer_start, crisprcrt_spacer_end, crisprdetect_dict[array]) # nearest crisprdetect row
							
							new_crisprdetect_row = construct_crisprdetect_row_minimal_update(row, nearest_row, crisprdetect_dict[array])
							new_array = new_crisprdetect_row[1]
							new_crisprdetect_row = new_crisprdetect_row[0]

							crisprdetect_dict[genome_id] = new_array
							crisprdetect_dict[genome_id].append(new_crisprdetect_row[0]) # add new row to show whether the row is CRT, CRISPRDETECT, or a joint prediction
							crisprdetect_dict_spacers[genome_id].add(new_crisprdetect_row[0][9]) # spacers

							break

						i += 1
		else:
		#	print("New CRISPR!!!")
			crisprdetect_dict[array] = add_array(crisprcrt_dict[array])
			# final code to add an entire new array goes here!!

	# add a union descriptor to any unmmodified rows.


	ret_out = open(output_dir, "w")
	spam_writer = csv.writer(ret_out)
	spam_writer.writerow(["Genome_id","orientation","orientation_score", "orientation_confidence", "questionable_array", "array_score","CRISPR-start","CRISPR-end", "repeat_start", "repeat_end","spacer_start","spacer_end","dr_repeat_original", "dr_repeat_concensous", "spacer", "Array_tool" ])
	outfiles = crisprdetect_dict.values() # this should create a list of lists
	for out_arr in outfiles:
		for row in out_arr:
			spam_writer.writerow(row)




	# write flattened_crisprdetect_dict.

			# add a new array in CRISPRdetect format!!	
			# need to develop a function for doing this

	# 1. Need to organise the entries by array and iterate.
	# 2. Compare array to array
	# 3. If the spacers from one array are in another array.
	# 	ignore.
	# Else if the spacers are not in the array but span the same region - ignore
	# Else if the spacres are upstream or downstream then add the new spacers using the concensous repeat of the existing CRISPR array
	# Else if the array doesn't exists then add the new array <- this will maximise sensitivity at the expense of accuracy.
	# The alternative being to only take concensous matched spacers.
	# Should keep a note of which spacers were detected by both CRISPR detection tools and which were unique. Create a new category!!
#	crisprdetect_table_url.close()
#	crisprcrt_table_url.close()
	ret_out.close()
	return 0
