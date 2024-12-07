# script to extract the pilercr id and start/end sites in a usable form for filtration

from Bio import SeqIO  
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord
import sys
import re
import csv



def line_denuller (split_sequences):
	ret_list = []
	for ele in split_sequences:
		if (ele != ''):
			ret_list.append(ele)
	return ret_list

def pos_extractor(argx):
	first_array = re.compile('Array 1')	# may need re.DOTALL enabled?
	no_crisprs = re.compile("0 putative CRISPR arrays found")
	characters = re.compile('[0-9]') # these patterns need to be improved to be more bombproof
	equals = re.compile('=')
	repeat_pattern = re.compile('')
	sequence_file = open(argx, "r") # split on >
	sequence_mega_str = sequence_file.read()
	if (re.search(first_array, sequence_mega_str)):
		sequence_mega_str = sequence_mega_str.split("Array 1")
		sequence_mega_str = "".join(sequence_mega_str[1:])
		sequence_mega_str = sequence_mega_str.split('SUMMARY BY SIMILARITY') # This should always exist!! This will ONLY work as intended if split is case sensitive!!!!!
		sequence_mega_str = sequence_mega_str[0]
	elif (re.search(no_crisprs, sequence_mega_str)):
		print("no_crisprs_detected!!")
		return 1
	else: 
		print("error input mega string unrecognised!! Possible bug?")
		return 2	

	#now should start with a prefiltered file.

	sequences = sequence_mega_str.split(">") # hopefully the strings include \n characters. Loading the whole str uses a fair amount of memory
	ret_file = open(argx + "_real_arr_positions.csv", "a")
	switch = False
	sequences = sequences[1:]
	for sequence in sequences:
		sequence_lines = sequence.split("\n")
		line = sequence_lines[0]
		sequence_lines = sequence_lines[1:]
		header = line # first line
		ret_lines = []
		for line in sequence_lines:
			line_items = line.split(" ")
			line_items = line_denuller(line_items) # get rid of the empty elements in the list from the split.
			ret_line = []
			if (line_items != []):
				pass
			if (line_items == []):
				continue
				
			elif (switch == False and re.search(equals , line_items[0][0])): # if line[0] contains an equals
				switch = True
				continue
			elif (switch == True and re.search(characters, line_items[0][0])):
				ret_line.append(line_items[0])
				ret_line.append(line_items[1])
				ret_lines.append(ret_line)
				continue
			elif (switch == True and re.search(equals, line_items[0][0])):
				switch = False
				break
			else:
				continue	
		start_pos = ret_lines[0][0]
		end_pos = str(int(ret_lines[-1][1]) + int(ret_lines[-1][0]))
		spam_writer = csv.writer(ret_file)
		# writing in form: [ genome_id, crispr_start_pos, crispr_end_p]
		spam_writer.writerow(["".join(header)] + [start_pos] + [end_pos])
	ret_file.close()
	return 0

def true_repeat(dr_seq, dr_dot, deletions=True):
	# direct repeats and dots should be the same length. Need to go through the dots and substitute characters if non-dot character is found
	i = 0
	output_seq = list(dr_seq)
	dr_dot = list(dr_dot)
	while (i < len(dr_dot)):

		if (dr_dot[i] != "."):
			if(dr_dot[i] == '-'): # deletion
				if (deletions):
					output_seq.pop(i)
					dr_dot.pop(i)
				
			else:
				output_seq[i] = dr_dot[i]
		i += 1
	return "".join(output_seq)				



def all_info_extractor_tony_DR_spacers(argx):
	ret_list = []
	first_array = re.compile('Array 1')	# may need re.DOTALL enabled?
	no_crisprs = re.compile("0 putative CRISPR arrays found")
	characters = re.compile('[0-9]') # these patterns need to be improved to be more bombproof
	equals = re.compile('=')
	repeat_pattern = re.compile('') # need the right pattern. Maybe Sungyeon/Jeremy could help?
	sequence_file = open(argx, "r") # split on >
	sequence_mega_str = sequence_file.read()
	if (re.search(first_array, sequence_mega_str)):
		sequence_mega_str = sequence_mega_str.split("Array 1")
		sequence_mega_str = "".join(sequence_mega_str[1:])
		sequence_mega_str = sequence_mega_str.split('SUMMARY BY SIMILARITY') # This should always exist!! This will ONLY work as intended if split is case sensitive!!!!!
		sequence_mega_str = sequence_mega_str[0]
	elif (re.search(no_crisprs, sequence_mega_str)):
		print("no_crisprs_detected!!")
		return 1
	else: 
		print("error input mega string unrecognised!! Possible bug?")
		return 2

	sequences = sequence_mega_str.split(">") # hopefully the strings include \n characters. Loading the whole str uses a fair amount of memory
	ret_file = open(argx + "_full_real_arr_positions_tony.csv", "w")
	spam_writer = csv.writer(ret_file)
	spam_writer.writerow(["Genome_id","CRISPR_start_position","CRISPR_end_position","Spacer_start","Spacer_end","Repeat_start","Repeat_end","Repeat_size","left_flanking_sequence","Direct_repeat","Direct_repeat_concensous","Spacer"])	
	switch = False
	#print(sequences)
	sequences = sequences[1:]
	for sequence in sequences:
	#	print(sequence)
		sequence_lines = sequence.split("\n")
		line = sequence_lines[0]
		sequence_lines = sequence_lines[1:]
		header = line # first line
		ret_lines = []
		#print(sequence_lines)
		i = 0
		for line in sequence_lines:
			line_items = line.split(" ")
			line_items = line_denuller(line_items) # get rid of the empty elements in the list from the split.
			ret_line = []
			if (line_items != []):
				pass
			if (line_items == []):
				i += 1
				continue
				
			elif (switch == False and re.search(equals , line_items[0][0])): # if line[0] contains an equals
				switch = True
				i += 1
				continue
			elif (switch == True and re.search(characters, line_items[0][0])):
				if (len(line_items) < 4):
					i += 1
					continue
				ret_line.append(line_items[0])
				ret_line.append(line_items[1])
				ret_line.append(line_items[3])
				ret_line.append(line_items[-3])
				ret_line.append(line_items[-2])
				ret_line.append(line_items[-1])
				ret_lines.append(ret_line)
				# here might need to reimagine this list.
				i += 1
				continue
			elif (switch == True and re.search(equals, line_items[0][0])):
				switch = False
				final_line = sequence_lines[-1]
				line_items = line.split(" ")
				line_items = line_denuller(line_items) # get rid of the empty elements in the list from the split.
				ret_line = []
				ret_line.append(line_items[3])
				ret_lines.append(ret_line)
				final_sequence_line = sequence_lines[i + 1]
				break
			else:
				i += 1
				continue	

		# now need to find the start and end positions for each line!!
		i = 0
		dr_coords = []
		spacer_coords = []
		line_coord = []
		final_line_items = line_denuller(final_sequence_line)
		final_line_items = "".join(i for i in final_line_items if not i.isdigit()) # really shouldn't use this for readability reasons. But it is cool!!
		final_line_items = final_line_items.strip('')
		final_line_items = final_line_items.strip()
	#	print(final_line_items)
		dr_sequence = final_line_items
		start_pos = ret_lines[0][0]
		if (ret_lines[-1][-1][0] != "=" or len(ret_lines) == 1): # special filtering case for corrupted data!!
			switch = False
			continue

		arg1 = ret_lines[-2][0]
		arg2 = ret_lines[-2][1]
		arg3 = ret_lines[-2][-1]
		end_pos = str(int(arg1) + int(arg2) + int(len(arg3)))
		while (i < len(ret_lines) - 1): #exclude final line
			spacer_start_coord = int(ret_lines[i][0]) + int(ret_lines[i][1]) # start_pos + DR
			spacer_end_coord = str(int(spacer_start_coord) + len(ret_lines[i][-1])) # length of spacer
			dr_start_coord = ret_lines[i][0]
			dr_end_coord = spacer_start_coord
			spacer_sequence = ret_lines[i][-1] # spacer sequence
			if ("." in spacer_sequence):
				i += 1
				continue
			dr_dots = ret_lines[i][-2]
			print(dr_sequence)
			my_dr_sequence = true_repeat(dr_sequence, dr_dots)
			left_flanking_sequence = ret_lines[i][-3] # check this in testing. May be more trouble than it's worth to include this information.
			dr_size = str(int(dr_end_coord) - int(dr_start_coord))
			line_coord.append(["".join(header)] + [start_pos, end_pos, spacer_start_coord, spacer_end_coord, dr_start_coord, dr_end_coord, dr_size, left_flanking_sequence, dr_sequence, my_dr_sequence, spacer_sequence])
			spam_writer.writerow(line_coord[-1])
			ret_list.append(line_coord[-1])
			i += 1
		# 2nd last line		
		# writing in form: [ genome_id, crispr_start_pos, crispr_end_p]	
	ret_file.close()
	return ret_list

def all_info_extractor (argx):
	ret_list = []
	first_array = re.compile('Array 1')	# may need re.DOTALL enabled?
	no_crisprs = re.compile("0 putative CRISPR arrays found")
	characters = re.compile('[0-9]') # these patterns need to be improved to be more bombproof
	equals = re.compile('=')
	repeat_pattern = re.compile('') # need the right pattern. Maybe Sungyeon/Jeremy could help?
	sequence_file = open(argx, "r") # split on >
	sequence_mega_str = sequence_file.read()
	if (re.search(first_array, sequence_mega_str)):
		sequence_mega_str = sequence_mega_str.split("Array 1")
		sequence_mega_str = "".join(sequence_mega_str[1:])
		sequence_mega_str = sequence_mega_str.split('SUMMARY BY SIMILARITY') # This should always exist!! This will ONLY work as intended if split is case sensitive!!!!!
		sequence_mega_str = sequence_mega_str[0]
	elif (re.search(no_crisprs, sequence_mega_str)):
		print("no_crisprs_detected!!")
		return 1
	else: 
		print("error input mega string unrecognised!! Possible bug?")
		return 2

	sequences = sequence_mega_str.split(">") # hopefully the strings include \n characters. Loading the whole str uses a fair amount of memory
	ret_file = open(argx + "_full_real_arr_positions.csv", "w")
	spam_writer = csv.writer(ret_file)
	spam_writer.writerow(["Genome_id","CRISPR_start_position","CRISPR_end_position","Spacer_start","Spacer_end","Repeat_start","Repeat_end","Repeat_size","left_flanking_sequence","Direct_repeat","Direct_repeat_concensous","Spacer"])
	switch = False
	sequences = sequences[1:]
	for sequence in sequences:
		sequence_lines = sequence.split("\n")
		line = sequence_lines[0]
		sequence_lines = sequence_lines[1:]
		header = line # first line
		ret_lines = []
		i = 0
		for line in sequence_lines:
			line_items = line.split(" ")
			line_items = line_denuller(line_items) # get rid of the empty elements in the list from the split.
			ret_line = []
			if (line_items != []):
				pass

			if (line_items == []):
				i += 1
				continue
				
			elif (switch == False and re.search(equals , line_items[0][0])):
				switch = True
				i += 1
				continue
			elif (switch == True and re.search(characters, line_items[0][0])):
				if (len(line_items) < 4):
					i += 1
					continue
				ret_line.append(line_items[0])
				ret_line.append(line_items[1])
				ret_line.append(line_items[3])
				ret_line.append(line_items[-3])
				ret_line.append(line_items[-1])
				ret_lines.append(ret_line)
				# here might need to reimagine this list.
				i += 1
				continue
			elif (switch == True and re.search(equals, line_items[0][0])):
				switch = False
				final_line = sequence_lines[-1]
				line_items = line.split(" ")
				line_items = line_denuller(line_items) # get rid of the empty elements in the list from the split.
				ret_line = []
				ret_line.append(line_items[3])
				ret_lines.append(ret_line)
				final_sequence_line = sequence_lines[i + 1]
				break
			else:
				i += 1
				continue	

		# now need to find the start and end positions for each line!!
		i = 0
		dr_coords = []
		spacer_coords = []
		line_coord = []
		final_line_items = line_denuller(final_sequence_line)
		final_line_items = "".join(i for i in final_line_items if not i.isdigit()) # really shouldn't use this for readability reasons. But it is cool!!
		final_line_items = final_line_items.strip('')
		final_line_items = final_line_items.strip()
		dr_sequence = final_line_items
		start_pos = ret_lines[0][0]
		if (ret_lines[-1][-1][0] != "=" or len(ret_lines) == 1): # special filtering case for corrupted data!!
			switch = False
			continue
		arg1 = ret_lines[-2][0]
		arg2 = ret_lines[-2][1]
		arg3 = ret_lines[-2][-1]
		end_pos = str(int(arg1) + int(arg2) + int(len(arg3)))
		while (i < len(ret_lines) - 2): #exclude final line and penultimate line
			spacer_start_coord = int(ret_lines[i][0]) + int(ret_lines[i][1]) # start_pos + DR
			spacer_end_coord = str(int(spacer_start_coord) + len(ret_lines[i][-1])) # length of spacer
			dr_start_coord = ret_lines[i][0]
			dr_end_coord = spacer_start_coord
			spacer_sequence = ret_lines[i][-1] # spacer sequence
			if ("." in spacer_sequence):
				i += 1
				continue
			dr_dots = ret_lines[i][-2]
			my_dr_sequence = true_repeat(dr_sequence, dr_dots)	
			left_flanking_sequence = ret_lines[i][-2] # check this in testing. May be more trouble than it's worth to include this information.
			dr_size = str(int(dr_end_coord) - int(dr_start_coord))
			line_coord.append(["".join(header)] + [start_pos, end_pos, spacer_start_coord, spacer_end_coord, dr_start_coord, dr_end_coord, dr_size, left_flanking_sequence, my_dr_sequence, dr_sequence, spacer_sequence]) # need to keep original and consensous dr_sequence!!
			spam_writer.writerow(line_coord[-1])
			ret_list.append(line_coord[-1])
			i += 1
		spacer_start_coord = int(ret_lines[i][0]) + int(ret_lines[i][1]) # start_pos + DR
		spacer_end_coord = str(int(spacer_start_coord) + len(ret_lines[i][-1])) # length of spacer
		dr_start_coord = ret_lines[i][0]
		dr_end_coord = spacer_start_coord
		spacer_sequence = ret_lines[i][-1] # spacer sequence
		if ("." in spacer_sequence):
			i += 1
			continue
		dr_dots = ret_lines[i][-2]
		my_dr_sequence = true_repeat(dr_sequence, dr_dots)	
		left_flanking_sequence = ret_lines[i][-2] # check this in testing. May be more trouble than it's worth to include this information.
		dr_size = str(int(dr_end_coord) - int(dr_start_coord))
		line_coord.append(["".join(header)] + [start_pos, end_pos, spacer_start_coord, spacer_end_coord, dr_start_coord, dr_end_coord, dr_size, left_flanking_sequence, my_dr_sequence, dr_sequence, ""]) # need to keep original and consensous dr_sequence!!
		spam_writer.writerow(line_coord[-1])
		ret_list.append(line_coord[-1])	
		# 2nd last line		
		# writing in form: [ genome_id, crispr_start_pos, crispr_end_p]	
	ret_file.close()
	return ret_list

def spacer_table_2_fasta_with_dr_offset (table_name, offset): # offset = number of nucleotides from the DR to include
	offset = int(offset)
	with open(table_name, "r") as csvfile:
		spacer_table = list(csv.reader(csvfile))
	spacer_table = spacer_table[1:]	
	ret_list = []	
	for spacer in spacer_table:
		new_spacer = Seq(spacer[-2][(-1 * offset):] + spacer[-1] + spacer[-2][: (offset - 1)])
		genome_base_position = spacer[0].split("::")
	#	print(genome_base_position)
		genome_base_position = genome_base_position[1]
		genome_base_position = genome_base_position.split(":")
		genome_base_position = int(genome_base_position[0])
		spacer_id = spacer[0]
		spacer_id = spacer[0].split(" ")
		spacer_id = spacer_id[0]
		spacer_id = spacer_id + "|" + "spacer_start_pos:" + spacer[3] + "|" + "spacer_end_pos:" + spacer[4] + "|" + "global_start_pos:" + str(int(spacer[3]) + genome_base_position) + "|" + "global_end_pos:" + str(int(spacer[4]) + genome_base_position)
		new_spacer_record = SeqRecord(new_spacer)
		new_spacer_record.id = spacer_id
		new_spacer_record.description = ""
		ret_list.append(new_spacer_record)
	SeqIO.write(ret_list, table_name + "_drs_" + str(offset) + "_spacers.fasta", "fasta")	
	return 0
# generate spacers from a table of CRISPR array predictions with the adjacent Direct repeat ends barcoded onto the spacers (without CRISPRdetect)
def merged_detectcrtpilercr_table_2_fasta_with_dr_offset_all_cases (table_name, offset):
	offset = int(offset)
	with open(table_name, "r") as csvfile:
		spacer_table = list(csv.reader(csvfile))
	spacer_table = spacer_table[1:]	
	ret_list = []	
	for spacer in spacer_table:
		# This doesn't work in cases of direct repeat permutation from the concensous
		new_spacer = Seq(spacer[12][(-1 * offset):] + spacer[14] + spacer[12][: (offset - 1)])
		new_spacer_left = Seq(spacer[12][(-1 * offset):] + spacer[14])
		new_spacer_right = Seq(spacer[14] + spacer[12][: (offset - 1)])
		genome_base_position = spacer[0].split("::")
		genome_base_position = genome_base_position[1]
		genome_base_position = genome_base_position.split(":")
		genome_base_position = int(genome_base_position[0])
		spacer_id = spacer[0]
		spacer_id = spacer[0].split(" ")
		spacer_id = spacer_id[0]
		spacer_id = spacer_id + "|" + "spacer_start_pos:" + spacer[10] + "|" + "spacer_end_pos:" + spacer[11] + "|" + "global_start_pos:" + str(int(spacer[10]) + genome_base_position) + "|" + "global_end_pos:" + str(int(spacer[11]) + genome_base_position)
		new_spacer_record = SeqRecord(new_spacer)
		new_spacer_record.id = spacer_id
		new_spacer_record.description = ""
		ret_list.append(new_spacer_record)

		new_spacer_record_left = SeqRecord(new_spacer_left)
		new_spacer_record_right = SeqRecord(new_spacer_right)
		new_spacer_record_left.id = spacer_id
		new_spacer_record_right.id = spacer_id 
		new_spacer_record_left.description = ""
		new_spacer_record_right.description = ""
		ret_list.append(new_spacer_record_left)
		ret_list.append(new_spacer_record_right)
	SeqIO.write(ret_list, table_name + "_drs_" + str(offset) + "_reconciled_spacers.fasta", "fasta")	
	return 0

# generate spacers from a table of CRISPR array predictions with the adjacent Direct repeat ends barcoded onto the spacers.
def merged_detectcrtpilercr_table_2_fasta_with_dr_offset_all_cases_but_better (table_name, offset):
	offset = int(offset)
	with open(table_name, "r") as csvfile:
		spacer_table = list(csv.reader(csvfile))

	ret_list = []	

	crisprs = {}
	for spacer in spacer_table:
		crispr = spacer[0] + spacer[1] + spacer[6] + spacer[7] # use split to take in the genome id and nothing more!!
		if (crispr not in crisprs):
			crisprs[crispr] = [spacer] # This method relies heavily on the crispr start/end being correct and the 
		else:
			crisprs[crispr].append(spacer)
	crisprs = crisprs.values()
	for crispr in crisprs:
		i = 0
		while (i < len(crispr) - 1):
			if (crispr[i][-2] == "" or crispr[i][10] == "" or crispr[i][11] == ""):
				i += 1	
				continue
			new_spacer = Seq(crispr[i][12][(-1 * offset):] + crispr[i][14] + crispr[i+1][12][: (offset - 1)])
			new_spacer_left = Seq(crispr[i][12][(-1 * offset):] + crispr[i][14])
			new_spacer_right = Seq(crispr[i][14] + crispr[i+1][12][: (offset - 1)])
			genome_base_position = crispr[i][0].split("::")
			genome_base_position = genome_base_position[1]
			genome_base_position = genome_base_position.split(":")
			genome_base_position = int(genome_base_position[0])
			
			spacer_id = crispr[i][0]
			spacer_id = crispr[i][0].split(" ")
			spacer_id = spacer_id[0]
			spacer_id = spacer_id + "|" + "spacer_start_pos:" + crispr[i][10] + "|" + "spacer_end_pos:" + crispr[i][11] + "|" + "global_start_pos:" + str(int(crispr[i][10]) + genome_base_position) + "|" + "global_end_pos:" + str(int(crispr[i][11]) + genome_base_position) + "|" + "array_tool:" + str(crispr[i][15])
			new_spacer_record = SeqRecord(new_spacer)
			new_spacer_record.id = spacer_id
			new_spacer_record.description = ""
			ret_list.append(new_spacer_record)
			new_spacer_record_left = SeqRecord(new_spacer_left)
			new_spacer_record_right = SeqRecord(new_spacer_right)
			new_spacer_record_left.id = spacer_id
			new_spacer_record_right.id = spacer_id 
			new_spacer_record_left.description = ""
			new_spacer_record_right.description = ""
			ret_list.append(new_spacer_record_left)
			ret_list.append(new_spacer_record_right)
			i += 1
		# The last spacer-dr sequence should be from the concensous repeat as the modified repeat is absent.
	SeqIO.write(ret_list, table_name + "_drs_" + str(offset) + "_reconciled_spacers.fasta", "fasta")	
	return 0
# generate spacers from a table of CRISPR array predictions with the adjacent Direct repeat ends barcoded onto the spacers (without CRISPRdetect)
def spacer_table_2_fasta_with_dr_offset_all_cases (table_name, offset): # offset = number of nucleotides from the DR to include
	offset = int(offset)
	with open(table_name, "r") as csvfile:
		spacer_table = list(csv.reader(csvfile))
	spacer_table = spacer_table[1:]	
	ret_list = []	
	for spacer in spacer_table:
		new_spacer = Seq(spacer[-2][(-1 * offset):] + spacer[-1] + spacer[-2][: (offset - 1)])
		new_spacer_left = Seq(spacer[-2][(-1 * offset):] + spacer[-1])
		new_spacer_right = Seq(spacer[-1] + spacer[-2][: (offset - 1)])
		genome_base_position = spacer[0].split("::")
		genome_base_position = genome_base_position[1]
		genome_base_position = genome_base_position.split(":")
		genome_base_position = int(genome_base_position[0])

		#the same id will be appended to all sequences!! Check this does not cause issues!!
		spacer_id = spacer[0]
		spacer_id = spacer[0].split(" ")
		spacer_id = spacer_id[0]
		spacer_id = spacer_id + "|" + "spacer_start_pos:" + spacer[3] + "|" + "spacer_end_pos:" + spacer[4] + "|" + "global_start_pos:" + str(int(spacer[3]) + genome_base_position) + "|" + "global_end_pos:" + str(int(spacer[4]) + genome_base_position)
		new_spacer_record = SeqRecord(new_spacer)
		new_spacer_record.id = spacer_id
		new_spacer_record.description = ""
		ret_list.append(new_spacer_record)

		new_spacer_record_left = SeqRecord(new_spacer_left)
		new_spacer_record_right = SeqRecord(new_spacer_right)
		new_spacer_record_left.id = spacer_id
		new_spacer_record_right.id = spacer_id 
		new_spacer_record_left.description = ""
		new_spacer_record_right.description = ""
		ret_list.append(new_spacer_record_left)
		ret_list.append(new_spacer_record_right)


	SeqIO.write(ret_list, table_name + "_drs_" + str(offset) + "_spacers.fasta", "fasta")	
	return 0

def spacer_table_2_fasta_with_dr_offset_left (table_name, offset): # offset = number of nucleotides from the DR to include
	offset = int(offset)
	with open(table_name, "r") as csvfile:
		spacer_table = list(csv.reader(csvfile))
	spacer_table = spacer_table[1:]	
	ret_list = []	
	for spacer in spacer_table:
		new_spacer = Seq(spacer[-2][(-1 * offset):] + spacer[-1])
		genome_base_position = spacer[0].split("::")
		genome_base_position = genome_base_position[1]
		genome_base_position = genome_base_position.split(":")
		genome_base_position = int(genome_base_position[0])
		spacer_id = spacer[0]
		spacer_id = spacer[0].split(" ")
		spacer_id = spacer_id[0]
		spacer_id = spacer_id + "|" + "spacer_start_pos:" + spacer[3] + "|" + "spacer_end_pos:" + spacer[4] + "|" + "global_start_pos:" + str(int(spacer[3]) + genome_base_position) + "|" + "global_end_pos:" + str(int(spacer[4]) + genome_base_position)
		new_spacer_record = SeqRecord(new_spacer)
		new_spacer_record.id = spacer_id
		new_spacer_record.description = ""
		ret_list.append(new_spacer_record)
	SeqIO.write(ret_list, table_name + "_drs_left_" + str(offset) + "_spacers.fasta", "fasta")	
	return 0	

def spacer_table_2_fasta_with_dr_offset_right (table_name, offset): # offset = number of nucleotides from the DR to include
	offset = int(offset)
	with open(table_name, "r") as csvfile:
		spacer_table = list(csv.reader(csvfile))
	spacer_table = spacer_table[1:]	
	ret_list = []	
	for spacer in spacer_table:
		new_spacer = Seq(spacer[-1] + spacer[-2][: (offset - 1)])
		genome_base_position = spacer[0].split("::")
	#	print(genome_base_position)
		genome_base_position = genome_base_position[1]
		genome_base_position = genome_base_position.split(":")
		genome_base_position = int(genome_base_position[0])
		spacer_id = spacer[0]
		spacer_id = spacer[0].split(" ")
		spacer_id = spacer_id[0]
		spacer_id = spacer_id + "|" + "spacer_start_pos:" + spacer[3] + "|" + "spacer_end_pos:" + spacer[4] + "|" + "global_start_pos:" + str(int(spacer[3]) + genome_base_position) + "|" + "global_end_pos:" + str(int(spacer[4]) + genome_base_position)
		new_spacer_record = SeqRecord(new_spacer)
		new_spacer_record.id = spacer_id
		new_spacer_record.description = ""
		ret_list.append(new_spacer_record)
	SeqIO.write(ret_list, table_name + "_drs_right_" + str(offset) + "_spacers.fasta", "fasta")	
	return 0

# generate a set of spacers in FASTA format from csv table containing compiled array predictions
def merged_detectcrtpilercr_table_2_fasta (table_name):
	with open(table_name, "r") as csvfile:
		spacer_table = list(csv.reader(csvfile))
#	spacer_table = spacer_table[1:]
	ret_list = []
	for spacer in spacer_table:
		if (spacer[-2] == '' or spacer[10] == '' or spacer[11] == ''):
			continue
		new_spacer = Seq(spacer[-2])
		genome_base_position = spacer[0].split("::")
		genome_base_position = genome_base_position[1]
		genome_base_position = genome_base_position.split(":")
		genome_base_position = int(genome_base_position[0])
		spacer_id = spacer[0]
		spacer_id = spacer[0].split(" ") # may not need this line anymore! Unlikely to break code though.
		spacer_id = spacer_id[0]

		if (spacer[10] == ''):
			print(spacer)
			print(spacer_id)
		spacer_id = spacer_id + "|" + "spacer_start_pos:" + spacer[10] + "|" + "spacer_end_pos:" + spacer[11] + "|" + "global_start_pos:" + str(int(spacer[10]) + genome_base_position) + "|" + "global_end_pos:" + str(int(spacer[11]) + genome_base_position) + "|" + "array_tool:" + str(spacer[15])
		new_spacer_record = SeqRecord(new_spacer)
		new_spacer_record.id = spacer_id
		new_spacer_record.description = ""
		ret_list.append(new_spacer_record)
	SeqIO.write(ret_list, table_name + "_reconciled_spacers.fasta", "fasta")
	return 0	

# generate a set of spacers in FASTA format from csv table containing compiled array predictions (without CRISPRdetect)
def spacer_table_2_fasta (table_name): # need to think about whether it would be more efficent to parse in the table directly? This could be a program option??
	with open(table_name, "r") as csvfile:
		spacer_table = list(csv.reader(csvfile))
	spacer_table = spacer_table[1:]	
	ret_list = []
	for spacer in spacer_table:
		new_spacer = Seq(spacer[-1])
		genome_base_position = spacer[0].split("::")
	#	print(genome_base_position)
		genome_base_position = genome_base_position[1]
		genome_base_position = genome_base_position.split(":")
		genome_base_position = int(genome_base_position[0])
		spacer_id = spacer[0]
		spacer_id = spacer[0].split(" ")
		spacer_id = spacer_id[0]
		spacer_id = spacer_id + "|" + "spacer_start_pos:" + spacer[3] + "|" + "spacer_end_pos:" + spacer[4] + "|" + "global_start_pos:" + str(int(spacer[3]) + genome_base_position) + "|" + "global_end_pos:" + str(int(spacer[4]) + genome_base_position)
		new_spacer_record = SeqRecord(new_spacer)
		new_spacer_record.id = spacer_id
		new_spacer_record.description = ""
		ret_list.append(new_spacer_record)
	SeqIO.write(ret_list, table_name + "_spacers.fasta", "fasta")	
	return 0