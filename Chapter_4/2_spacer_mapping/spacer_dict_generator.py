# spacer_dict_generator

def spacer_dict_gen(argx):
	with open(argx) as csvfile2:
		spacer_table = list(csv.reader(csvfile2))

	spacer_table.sort(key=first_ele) # check this works in isolation
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
		spacer_start_pos = spacer_hit[0].split("spacer_start_pos:") # may be better to split on global_start/end_pos instead!!!
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
				coord_dict[spacer_coord] = [spacer_hit] # may be creating a superflous nested list of one element. Easy to bypass however!!
			else:
				coord_dict[spacer_coord].append(spacer_hit)
		else:
			spacer_dict[previous_spacer_id] = coord_dict
			coord_dict = {}
		previous_spacer_id = spacer_id	
	spacer_dict[previous_spacer_id] = coord_dict # should this be in the outer scope (as it was originally?) Possible error?
	return spacer_dict
