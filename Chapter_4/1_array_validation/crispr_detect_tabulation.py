# CRISPR_detect_tabulation

# program to tabulate CRISPRdetect results into csv

import csv
import sys
from Bio import SeqIO
import re
import pilercr_pos_extractor_annotation_full_spacers
import itertools

# function to adjust the tabulated positions and sequences of CRISPR direct repeats containing INDELS
def true_repeat_indels (repeat, dots, repeat_start, repeat_size, sense, indels_positions=[]):
	ori_repeat = list(repeat)
	ori_dots = list(dots)
	if (indels_positions == []):
		return [pilercr_pos_extractor_annotation_full_spacers.true_repeat(repeat,dots),0]
	else:
		if (re.search("Deletion", "".join(indels_positions))):
			return (8,8)
		else:	

			indel_info = indels_positions	
			indels = indels_positions
			mutations = indels[0]			
			coordinates  = indels [1]
			mutations = mutations.split(",")
			coordinates = coordinates[1:-1].split(",")
			repeat_vector = []
			if sense == 1:
				for coord in coordinates:
					repeat_coord = int(coord) - int(repeat_start) 
					repeat_vector.append(repeat_coord)
			else:
				for coord in coordinates:
					repeat_coord = int(repeat_start) - int(coord)
					repeat_vector.append(repeat_coord)
			i = 0
			while (i < len(mutations)): # should be the same size as the repeat vector
			#	 print(mutations[i], repeat_vector[i])
				 ori_repeat.insert(repeat_vector[i], mutations[i])
				 ori_dots.insert(repeat_vector[i], list(len(mutations[i]) * "."))
				 i += 1

			flat_ori_repeat = list(itertools.chain(*ori_repeat))
			flat_ori_dots = list(itertools.chain(*ori_dots))			
			mutation_str = "".join(mutations)

			return [pilercr_pos_extractor_annotation_full_spacers.true_repeat(flat_ori_repeat,flat_ori_dots), len(mutation_str)] 



# function to convert the output file produced by CRISPRdetect into csv format.
def tabulate(crispr_detect_url):
	
	ret_out = open(crispr_detect_url + "_crisprdetect_results.csv", "w")
	spam_writer = csv.writer(ret_out)
	spam_writer.writerow(["genome_id","orientation","orientation_score", "orientation_confidence", "questionable_array", "array_score","CRISPR-start","CRISPR-end", "repeat_start", "repeat_end","spacer_start","spacer_end","dr_repeat_original", "dr_repeat_concensous", "spacer" ]) # add a specific crispr_id? Need to consider whether or not to adjust for CRISPR-position
	# INPUT: CRISPR-detect output file
	crisprdetect_handle = open(crispr_detect_url, "r")
	crisprdetect_results = crisprdetect_handle.read()
	crisprdetect_results = crisprdetect_results.split(">")
	crisprdetect_handle.close()
	for crispr in crisprdetect_results[1:]:
		ret_crispr = []
		my_crispr = crispr.split("\n")
		my_crispr = pilercr_pos_extractor_annotation_full_spacers.line_denuller(my_crispr)
		crispr_header = my_crispr[0]
		crispr_identifer = pilercr_pos_extractor_annotation_full_spacers.line_denuller (crispr_header.split("\t"))
		crispr_orientation = crispr_identifer[1]
		crispr_orientation.strip()
		crispr_orientation = crispr_orientation.split("Array_Orientation: ") [1]

		crispr_identifer = crispr_identifer[0]
		i = 1
		data_switch = 0
		dr_repeat = ""
		# iterate through the rows in the crispr
		while (i < len(my_crispr)):
			my_row = my_crispr[i].split()			
			my_row = pilercr_pos_extractor_annotation_full_spacers.line_denuller(my_row) # remove empty entries

			if (my_row[0] == "Position"):
				i += 2
				data_switch = 1
				continue
			if (my_row[0][0] == '=' and data_switch == 1):
				data_switch = 0
				dr_repeat = "".join(pilercr_pos_extractor_annotation_full_spacers.line_denuller( my_crispr[i + 1].split() [4] ))
				questionable_array = "".join(pilercr_pos_extractor_annotation_full_spacers.line_denuller( my_crispr[i + 4].split() [4] ))
				array_score = "".join(pilercr_pos_extractor_annotation_full_spacers.line_denuller( my_crispr [i + 4].split() [6] ))

				final_direction_score = "".join(pilercr_pos_extractor_annotation_full_spacers.line_denuller( my_crispr[i + 18].split() [4] ))
				final_direction_score = final_direction_score[1:].split(",")[0]
				final_direction_confidence = "".join(pilercr_pos_extractor_annotation_full_spacers.line_denuller( my_crispr[i + 18].split() [6]) [:-1])
				
				break
			if (data_switch == 1):
				ret_crispr.append(my_row)
			i += 1	
		k = 0 
		while (k < len(ret_crispr)):
			sense = 1
			if (crispr_orientation == "Reverse"):
				sense = -1

			crispr_start = ret_crispr[0][0]
			crispr_end = str(int(ret_crispr [-1][0]) + int(ret_crispr[-1][1]) * sense)
			repeat_start = ret_crispr[k][0]
			snv_dr = pilercr_pos_extractor_annotation_full_spacers.true_repeat(dr_repeat, ret_crispr[k][4],False)
			mutation_length = 0
			if (len(ret_crispr[k]) > 6 or '-' in ret_crispr[k][4]):
				if (len(ret_crispr[k]) <= 6):
					my_out = true_repeat_indels(snv_dr,ret_crispr[k][4], repeat_start, ret_crispr[k][3],sense)
					true_dr = my_out[0]
					mutation_length = my_out[1]
				else:
					my_out = true_repeat_indels(snv_dr,ret_crispr[k][4], repeat_start, ret_crispr[k][3],sense, ret_crispr[k][6:])
					true_dr = my_out[0]
					mutation_length = my_out[1]
			else:
				true_dr = snv_dr
	
			repeat_end = str(int(repeat_start) + (int(ret_crispr[k][1]) + mutation_length )* sense)

			spacer_start = str(int(repeat_end) + 1 * sense)
			spacer_end = str(int(spacer_start) + (int(ret_crispr[k][3]) - 1) * sense)

			# now just need to reconstitute the original drs.
			
			spacer = ret_crispr[k][5]

			if (true_dr == 8): # truncated spacer!! Consider whether to include this!!
				k += 1
				continue

			if ("-" in spacer): 
				k += 1
				continue	
			spam_writer.writerow([crispr_identifer, crispr_orientation, final_direction_score, final_direction_confidence,questionable_array,array_score,crispr_start,crispr_end, repeat_start, repeat_end,spacer_start, spacer_end,true_dr,dr_repeat,spacer ])

			k += 1

	ret_out.close()
	return 0	
	