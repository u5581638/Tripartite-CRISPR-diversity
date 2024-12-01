# crisprdetect_inversion

# function to invert crisprdetect output table so that all entry are sorted smallest to largest by gene numbers

import sys
import csv
import copy
from collections import OrderedDict

# helper function for sorting
def ret_value(x):
	return int(x[8])

# function to find reverse complement of a sequence
def rev_comp(sequence):
	ret_str = []
	for sequ in sequence:
		if (sequ == "a" or sequ == "A"):
			ret_str.append("T")
		elif (sequ == "t" or sequ == "T"):
			ret_str.append("A")	
		elif (sequ == "c" or sequ == "C"):
			ret_str.append("G")
		elif (sequ == "G" or sequ == "g"):
			ret_str.append("C")
		elif (sequ == "U" or sequ == "u"):
			ret_str.append("A")
		else:
			ret_str.append(sequ)
	ret_str.reverse()		
	return "".join(ret_str)							

# main function to invert the coordinates, spacers and Direct repeat sequences for arrays predicted on the antisense strand.
def invert(crisprdetect_table_url):
	# input table of CRISPR-detect predictions in csv format.
	crisprdetect_table_handle = open(crisprdetect_table_url, "r")
	crisprdetect_table = list(csv.reader(crisprdetect_table_handle))
#	print(crisprdetect_table)
	ret_out = open(crisprdetect_table_url + "_inverted.csv", "w")
	spam_writer = csv.writer(ret_out)
	ret_dict = OrderedDict()
	i = 0
	for row in crisprdetect_table:
		if (i == 0):
			i = 1
			continue
		
		my_row = copy.deepcopy(row)
		genome_id = my_row[0].split(" ") [0] # genome_id only

		if (my_row[1] == "Reverse"):
			spacer_start_cpy = copy.deepcopy(my_row[11])
			spacer_end_cpy = copy.deepcopy(my_row[10])
			repeat_start_cpy = copy.deepcopy(my_row[9])
			repeat_end_cpy = copy.deepcopy(my_row[8])
			array_end_cpy = copy.deepcopy(my_row[6])
			arr_start_cpy = copy.deepcopy(my_row[7])
			my_row[6] = arr_start_cpy
			my_row[7] = array_end_cpy
			my_row[10] = spacer_start_cpy
			my_row[11] = spacer_end_cpy
			my_row[8] = repeat_start_cpy
			my_row[9] = repeat_end_cpy
			my_row[12] = rev_comp(copy.deepcopy(my_row[12]))
			my_row[13] = rev_comp(copy.deepcopy(my_row[13]))
			my_row[14] = rev_comp(copy.deepcopy(my_row[14]))
		if (genome_id not in ret_dict):
			ret_dict[genome_id] = [my_row]
		else:	
			ret_dict[genome_id].append(my_row)	
#	ret_dict = list(ret_dict.values()).sort(key= ret_value) # hopefully this is still in the order the keys were added!!
	ret_dict = list(ret_dict.values())
	for array in ret_dict:
		if (array[0][1] == "Reverse"):
			array.sort(key=ret_value)
			sorted_array = copy.deepcopy(array)
			i = 1
			previous_row_repeat_start = copy.deepcopy(sorted_array[0][8])
			previous_row_repeat_end = str(int(copy.deepcopy(sorted_array[0][9])))
			previous_row_repeat_con = copy.deepcopy(sorted_array[0][12])
			previous_row_repeat = copy.deepcopy(sorted_array[0][13])
			sorted_array[0][10] = str(int(sorted_array[0][10]))
			sorted_array[0][11] = str(int(sorted_array[0][11]))

			# repeats may need to be incremented to adjust for row.
			while (i < len(sorted_array)):
			#	print(sorted_array[i-1][8])
				previous_row_repeat_start_cpy = copy.deepcopy(previous_row_repeat_start)
				previous_row_repeat_end_cpy = copy.deepcopy(previous_row_repeat_end)
				previous_row_repeat_con_cpy = copy.deepcopy(previous_row_repeat_con)
				previous_row_repeat_cpy = copy.deepcopy(previous_row_repeat)

				previous_row_repeat_start = copy.deepcopy(sorted_array[i][8])
				previous_row_repeat_end = copy.deepcopy(sorted_array[i][9])
				previous_row_repeat_con = copy.deepcopy(sorted_array[i][12])
				previous_row_repeat = copy.deepcopy(sorted_array[i][13])

				sorted_array[i][8] = previous_row_repeat_start_cpy
				sorted_array[i][9] = previous_row_repeat_end_cpy
				sorted_array[i][12] = previous_row_repeat_con_cpy
				sorted_array[i][13] = previous_row_repeat_cpy

				sorted_array[i][10] = str(int(sorted_array[i][10]) + 1)
			#	if (i == 2):
			#		sorted_array[i][11] = str(int(copy.deepcopy(sorted_array[i][11])) + 1) 
			#	else:	
			#		sorted_array[i][11] = str(int(sorted_array[i][11]) + 1) 
				i += 1
			final_row = [sorted_array[i-1][0], sorted_array[i-1][1],sorted_array[i-1][2],sorted_array[i-1][3],sorted_array[i-1][4],sorted_array[i-1][5],sorted_array[i-1][6],sorted_array[i-1][7],previous_row_repeat_start,previous_row_repeat_end,"","",previous_row_repeat_con,previous_row_repeat,""]
			sorted_array = sorted_array[1:]
			sorted_array.append(final_row)
		else:
			sorted_array = array
			sorted_array[-1][10] = ""
			sorted_array[-1][11] = ""

		
		for row in sorted_array:
			if (row[14] == '|'):
				row[14] = "" # This should now be partially redundant
			spam_writer.writerow(row)	
	ret_out.close()
	crisprdetect_table_handle.close()
	return 0	
		





