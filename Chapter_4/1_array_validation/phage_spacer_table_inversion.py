import csv 
import sys
import pandas

def phage_table_invert(hitmap_table_url):

	with open(hitmap_table_url) as csvfile:
		hit_table = pandas.read_csv(csvfile)
	ret_dict = {}
	ret_dict['Phage_mapped_site'] = []
	ret_dict['Phage'] = []
	ret_dict['mapped_start_coordinate'] = []
	ret_dict['mapped_end_coordinate'] = []
	ret_dict['mapped_evalue'] = []
	ret_dict['genome'] = []
	ret_dict['genome_spacer_start'] = []
	ret_dict['genome_spacer_end'] = []
	ret_dict['genome_sense'] = []
	ret_dict['genome_global_start_pos'] = []
	ret_dict['genome_global_end_pos'] = []
	ret_dict['RUN'] = []
	for hit in hit_table:
		phage_id = hit ['mapped_phage_id'] # replace this with the column name for the mapped phage!! 
		mapped_start_coordinate = hit ['start_coord'] # BLAST identified start
		mapped_end_coordinate = hit ['end_coord']
	#	mapped_sense = should the sense implictly by the order of the start and end coordinates
		spacer_mapped_site = phage_id + "_" + mapped_start_coordinate + "_" + mapped_end_coordinate  # this should be the phage id + start_mapped_position + end_mapped_position 
		mapped_evalue = hit ['evalue']
		spacer_id = hit['spacer']
		spacer_id.split("|")
		genome = spacer_id[0]
		genome_spacer_start = spacer_id[1].split(":") [1]
		genome_spacer_end = spacer_id[2].split(":") [1]
		genome_sense = hit['orientation']
		global_start_pos = spacer_id[3].split(":")[1]
		global_end_pos = spacer_id[4].split(":") [1]
		run = hit['RUN']
		ret_dict['RUN'].append(run)
		ret_dict['Phage_mapped_site'].append( spacer_mapped_site)
		ret_dict['Phage'].append( phage_id)
		ret_dict['mapped_start_coordinate'].append( mapped_start_coordinate)
		ret_dict['mapped_end_coordinate'].append( mapped_end_coordinate)
		ret_dict['mapped_evalue'].append( mapped_evalue)
		ret_dict['genome'].append( genome )
		ret_dict['genome_spacer_start'].append( genome_spacer_start)
		ret_dict['genome_spacer_end'].append( genome_spacer_end)
		ret_dict['genome_sense'].append( genome_sense )
		ret_dict['genome_global_start_pos'].append( global_start_pos)
		ret_dict['genome_global_end_pos'].append(  global_end_pos) 


	return pandas.DataFrame(ret_dict)



