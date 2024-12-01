# pilercr_reconcilation
# reconcile pilercr_predictions with a table of crt/CRISPRdetect merged results

import sys
import csv
import copy


# BEST TO DO THIS TOMORROW WHEN WIDE AWAKE!!!!!!!
# I Suspect the final reconciliation of CRISPRs will be an all-by-all effort on Friday
# Failing this then over the weekend!!!
# Should be done by Monday!!!!!
# Need
def add_array (pilercr_array):
	ret_list = []
	dr_con = []
	for array in pilercr_array:
		ret_list.append([array[0].split(" ") [0], "Unknown","NA", "NA","NA","NA",array[1],array[2], array[5],array[6], array[3], array[4], array[9],array[10],array[11], "PILERCR"] )
		dr_con.append(array[10])
	con_dr = max(set(dr_con), key=dr_con.count)
	i = 0
	while (i < len(ret_list)):
		ret_list[i][13] = con_dr
		i += 1
	return ret_list

def nearest_array (pilercr_spacer_start, pilercr_spacer_end, array):	
	distances = {}
	for row in array:
		crispr_start = int(row[6])
		crispr_end = int(row[7])
		distance = (abs((pilercr_spacer_start - crispr_start) + (pilercr_spacer_end - crispr_end))) / 2
		distances[distance] = row
	shortest_distance = min(distances.keys())
	return distances[shortest_distance]	

# construct a new row with maximal direct repeat corrdinates to be consistent with predictions from CRISPR-CRT an CRISPRdetect
def construct_crisprdetect_row_minimal(row, nearest_row):
	ret_row = []
	spacer_start = int(row[5])
	spacer_end = int(row[6])
	start_distance = int(spacer_start) - int(nearest_row[6])
	end_distance = int(spacer_end) - int(nearest_row[7])

	if (row[1] < nearest_row[6]):
		crispr_start = row[1]
	else:
		crispr_start = nearest_row[6]
	if (row[2] > nearest_row[7]):
		crispr_end = row[2]
	else:
		crispr_end = nearest_row[7]
	ret_row = [[row[0].split(" ") [0]] + nearest_row[1:6] + [crispr_start] + [crispr_end] + [row[5]] + [row[6]] + [row[3]] + [row[4]] + [row[9]] + [nearest_row[13]] + [row[11]] + ["PILERCR"]]	
	return ret_row

# function to reconile piler-cr array predictions with those from reconciled CRISPRdetect/CRISPR-CRT
def pilercr_reconcile (crisprdetect_table_url, crispr_pilercr_url, output_dir):

	csvfile = open(crisprdetect_table_url, "r")
	csvfile2 = open(crispr_pilercr_url, "r")

	crisprdetect_table = list(csv.reader(csvfile))
	crispr_pilercr_table = list(csv.reader(csvfile2))

	csvfile.close()
	csvfile2.close()
	# 1. Dictionalise all entries into tabular form!!
	# CRISPR_detect table should be pre-inverted!!

	crisprdetect_dict = {}
	crisprdetect_dict_spacers = {}

	for row in crisprdetect_table[1:]:
		my_row = copy.deepcopy(row)
		genome_id = my_row[0].split(" ") [0] # genome_id only

		if genome_id not in crisprdetect_dict:
			crisprdetect_dict[genome_id] = [my_row]
			crisprdetect_dict_spacers[genome_id] = {my_row[14]}
		else:
			crisprdetect_dict[genome_id].append(my_row)
			crisprdetect_dict_spacers[genome_id].add(my_row[14])

	crispr_pilercr_dict = {}

	for row in crispr_pilercr_table[1:]:
		genome_id = row[0].split(" ") [0]
		if (genome_id not in crispr_pilercr_dict):
			crispr_pilercr_dict[genome_id] = [row] # should this be a dictionary indexes spacers to theur respective rows? No need to create a seperate dictionary doing this mapping and assign spacers as the values of 

		else:
			crispr_pilercr_dict[genome_id].append(row)

	for array in crispr_pilercr_dict:
		if (array in crisprdetect_dict):
			for row in crispr_pilercr_dict[array]:		
				if (row[11] in crisprdetect_dict_spacers[array]):
					i = 0
					while i < len (crisprdetect_dict[array]):
						# this needs to be dictionalised to prevent additional iteration!!!!!!!!!!!!!!!!!
						# need to be able to see all rows as list/set in order to be able to decide where to add new spacers. 
						# need a special condition if insertion results in an overlap with another CRISPR-array!!!!
						if (row[11] == crisprdetect_dict[array][i][14]):
							print("Key row!!")
							print(row[11])
							if (crisprdetect_dict[array][i][15] == "CRISPRDETECT"):
								crisprdetect_dict[array][i][15] = "CRISPRDETECT.PILERCR" # Add a note that spacers match to both
							elif (crisprdetect_dict[array][i][15] == "CRT"):
								crisprdetect_dict[array][i][15] = "CRT.PILERCR"
							elif (crisprdetect_dict[array][i][15] == "CRISPRDETECT.CRT"):
								crisprdetect_dict[array][i][15] = "CRISPRDETECT.CRT.PILERCR"
							else:
								print("Unrecognised entry:")
								print(crisprdetect_dict[array][i])	
							
						i += 1	 
				else:
					crispr_pilercr_spacer_start = int(row[1])
					crispr_pilercr_spacer_end = int(row[2])
					i = 0
					while i < len (crisprdetect_dict[array]):
						# this needs to be dictionalised to prevent additional iteration!!!!!!!!!!!!!!!!!
						# need to be able to see all rows as list/set in order to be able to decide where to add new spacers. 
						# need a special condition if insertion results in an overlap with another CRISPR-array!!!!
						crisprdetect_start = int(crisprdetect_dict[array][i][6])
						crisprdetect_end = int(crisprdetect_dict[array][i][7])
						genome_id = crisprdetect_dict[array][i][0].split(" ") [0]
						if (crisprdetect_start <= crispr_pilercr_spacer_start <= crisprdetect_end or crisprdetect_start <= crispr_pilercr_spacer_end <= crisprdetect_end): # if the spacer is within the array existing array
							# assume this is a pre-existing spacer.
							if (crisprdetect_dict[array][i][15] == "CRISPRDETECT"):
								crisprdetect_dict[array][i][15] = "CRISPRDETECT.PILERCR" # Add a note that spacers match to both
							elif (crisprdetect_dict[array][i][15] == "CRT"):
								crisprdetect_dict[array][i][15] = "CRT.PILERCR"
							elif (crisprdetect_dict[array][i][15] == "CRISPRDETECT.CRT"):
								crisprdetect_dict[array][i][15] = "CRISPRDETECT.CRT.PILERCR"
							else:
								print("Unrecognised entry:")
								print(crisprdetect_dict[array][i])	
							#	exit()
						else:
							# need to add to nearest array
							print(i)
							print(array)
							print(crisprdetect_start)
							print(crisprdetect_end)
						#	print(crisprdetect_dict[array])
							nearest_row = nearest_array(crispr_pilercr_spacer_start, crispr_pilercr_spacer_end, crisprdetect_dict[array]) # nearest crisprdetect row
							new_crisprdetect_row = construct_crisprdetect_row_minimal(row, nearest_row)
							print("Yay!!")
							print(new_crisprdetect_row)
							new_crisprdetect_row[0][15] = "PILERCR"
							crisprdetect_dict[genome_id].append(new_crisprdetect_row[0]) # add new row to show whether the row is CRT, CRISPRDETECT, or a joint prediction
						#	print(crisprdetect_dict_spacers[genome_id])
							crisprdetect_dict_spacers[genome_id].add(new_crisprdetect_row[0][14])

							break

						i += 1	

		else:
			crisprdetect_dict[array] = add_array(crispr_pilercr_dict[array])
			# final code to add an entire new array goes here!!

	# THINK ABOUT WHETHER THIS LOOP IS NESSESSARY, WHAT MISSING VALUES IT ASSIGNS

	ret_out = open(output_dir, "w")
	spam_writer = csv.writer(ret_out)
#	spam_writer.writerow(["genome_id","orientation","orientation_score", "orientation_confidence", "questionable_array", "array_score","CRISPR-start","CRISPR-end", "repeat_start", "repeat_end","spacer_start","spacer_end","dr_repeat_original", "dr_repeat_concensous", "spacer", "Array_tool" ])
	outfiles = crisprdetect_dict.values() # this should create a list of lists
	print("Ready to write!!")
	print(len(outfiles))
	for out_arr in outfiles:
		for row in out_arr:
			spam_writer.writerow(row)
	return 0				