# spacer_dict_generator

# generate a dictionary of spacers and their associated DRs with coordinates
def spacer_dict_gen(argx):
	with open(argx) as csvfile2:
		spacer_table = list(csv.reader(csvfile2))

	spacer_table.sort(key=first_ele) 
	spacer_dict = {}
	# want to create a two level nested dict
	# the outer level will be the genome id
	# the inner level will be the spacer start and end positions
	# the values will be the merger of the spacer_hits and spacer_info
	coord_dict = {}
	previous_spacer_id = ""

	# result dict for lookup of 
	for spacer_hit in spacer_table:
		# extract the genome_id
		spacer_id = spacer_hit[0].split("|")
		spacer_id = spacer_id[1]		
		spacer_start_pos = spacer_hit[0].split("spacer_start_pos:") 
		spacer_start_pos = spacer_start_pos[1]
		spacer_start_pos = spacer_start_pos.split("|")
		spacer_start_pos = spacer_start_pos[0]
		spacer_end_pos = spacer_hit[0].split("spacer_end_pos:")
		spacer_end_pos = spacer_end_pos[1]
		spacer_end_pos = spacer_end_pos.split("|")
		spacer_end_pos = spacer_end_pos[0]
		spacer_coord = (spacer_start_pos, spacer_end_pos)
	#	print(spacer_coord)
		if (spacer_id == previous_spacer_id):
			if (spacer_coord not in coord_dict):
				coord_dict[spacer_coord] = [spacer_hit] 
			else:
				coord_dict[spacer_coord].append(spacer_hit)
		else:
			spacer_dict[previous_spacer_id] = coord_dict
			coord_dict = {}
		previous_spacer_id = spacer_id	
	spacer_dict[previous_spacer_id] = coord_dict
	return spacer_dict
