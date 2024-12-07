from Bio import SeqIO
import sys
import csv
import re

# functions to parse CRISPRcrt/PILER-CR predictions into tables in csv format.

# helper function to dictionalise the strand direction info table produced by CRISPRleader
def dictionary_generator (info_table):
	crispr_info_dict = {}
	info_switch = 0
	for row in info_table:
		if (re.match('=', row[0][0]) and info_switch == 0):
			info_switch = 1
			continue
		if (re.match('=', row[0][0]) and info_switch == 1):
			info_switch = 0	
		if (info_switch == 1 and re.match('[0-9]', row[0][0])):
			my_row = row[0].split('\t')
			crispr_info_dict[my_row[0]] = my_row # should already be in list format		
	return crispr_info_dict
	# info dissection goes here!!

# function to tabulate output from CRISPR-CRT (as part of CRISPRleader) in csv format.
def spacer_crt_table_generation (spacers_table, info_table, output_dir):
	spacer_crt_table = {}
	crispr_info_dict = dictionary_generator(info_table)
	i=0
	switch = 0
	numbering = 0
	crispr_no_iter = 0
	while (i < len(spacers_table)):
		if (spacers_table[i] != []):
		#	print(spacers_table[i][0])
			if re.match("ORGANISM", spacers_table[i][0]):

				organism_line = spacers_table[i]
				organism_line = organism_line[0].split("ORGANISM:  ")
				organism_line = organism_line[1]
				organism_line = organism_line.split(" ")
				organism_line = organism_line[0]
			if re.match("CRISPR", spacers_table[i][0]):
				crispr_no_iter += 1
				crispr_no = str(crispr_no_iter)

			if (re.match('-', spacers_table[i][0][0]) and switch == 0):
				switch = 1
				my_crispr = []
				i += 1 
				continue

			if (re.match('-', spacers_table[i][0][0]) and switch == 1):
				switch = 0
				if(numbering == 1):
					k=0
					while (k < len(my_crispr)):
						my_crispr[k].append(str(k))
						k += 1
				elif(numbering == 2):
					k= len(my_crispr) - 1
					c = 0
					while (k >= 0):
						my_crispr[k].append(str(c))
						c += 1
						k -= 1
				else:
					k=0
					while (k < len(my_crispr)):
						my_crispr[k].append("N/A")
						k += 1
			
				if (my_crispr[0][0] in spacer_crt_table):
					spacer_crt_table[my_crispr[0][0]].extend(my_crispr) # Might need id simplification. Might want to check what this does!!
				else:
					spacer_crt_table[my_crispr[0][0]] = my_crispr	
					
			if (switch == 1 and re.match('[0-9]', spacers_table[i][0][0]) and len(spacers_table[i][0].split('\t')) >= 2): # originally 5 

				my_row = spacers_table[i][0].split("\t")
				repeat_start = my_row[0]
				if (len(spacers_table[i][0].split('\t')) >= 5):
					repeat_end = my_row[-1][1:]
				else:
					repeat_end = str(len(my_row[2]) -1)
				repeat_end_coord = repeat_end.split(",")
				repeat_end_coord_r = repeat_end_coord[0].split(" ")
				repeat_end_coord_r = repeat_end_coord_r[-1]
				repeat_end_coord_s = repeat_end_coord[0].split(" ")
				if (len(spacers_table[i][0].split('\t')) >= 5):
					repeat_end_coord_s = repeat_end_coord_s[1]
					spacer_start = int(repeat_start) + int(repeat_end_coord_r)
					spacer_end = spacer_start + len(my_row[3])
					row_line = [organism_line, "CRISPR_" + str(crispr_no), repeat_start, str(int(repeat_start) + int(repeat_end)), spacer_start, str(int(spacer_end)),my_row[-3], my_row[-2]]
				else:
					repeat_end_coord_s = ""
					spacer_start = ""
					spacer_end = ""
					row_line = [organism_line, "CRISPR_" + str(crispr_no), repeat_start, str(int(repeat_start) + int(repeat_end)), spacer_start, str(spacer_end),my_row[2], ""]
				if (crispr_no in crispr_info_dict):
					row_line.extend(crispr_info_dict[crispr_no])
					if (crispr_info_dict[crispr_no][-4] == "Forward"):
						numbering = 1
					elif (crispr_info_dict[crispr_no][-4] == "Reverse"):
						numbering = 2
					else:
						numbering = 0		

				my_crispr.append(row_line)

		i += 1	
	with open(output_dir + "_orientation_out.csv", "w") as csvfile2:
		spam_writer = csv.writer(csvfile2)
		spam_writer.writerow(["Genome_id","CRISPR_number","CRISPR-start_placeholder","CRISPR-end_placeholder","Repeat_start","Repeat_end","Spacer_start","Spacer_end","Direct_repeat","Spacer", "CRISPR_Array_number", "Direct_Repeat","CRISPR_array_range","Array_orientation", "Repeat_cluster","Leader_group","CRISPRleader","CRISPRleader_length"])
		my_table = list(spacer_crt_table.values())
		for row in my_table:
			for my_row in row:
				custom_row = my_row[10].split("-")
				array_start = custom_row[0]
				array_end = custom_row[1] 
				if (my_row[4] == ''):
					spam_writer.writerow(my_row[:2] + [array_start] + [array_end] + [str(int(my_row[2]) - 0)] + [my_row[3]] + [""] + [""] + my_row[6:])
				else:
					ret_row = my_row[:2] + [array_start] + [array_end] + [str(int(my_row[2]) - 0)] + [my_row[3]] + [str(int(my_row[4]) + 1)] + [str(int(my_row[5]))] + my_row[6:]
					spam_writer.writerow(my_row[:2] + [array_start] + [array_end] + [str(int(my_row[2]) - 0)] + [my_row[3]] + [str(int(my_row[4]) + 1)] + [str(int(my_row[5]))] + my_row[6:])
	return spacer_crt_table
