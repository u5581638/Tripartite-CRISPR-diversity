# spacers 2 or more

import csv
import sys

# unused function. 
def filter(spacer_url, output_url):
	with open (spacer_url, "r") as csvfile:
		hit_table = list(csv.reader(csvfile))
	spacer_mapped_dict = {}
	outf = open(output_url,"w")
#	header = hit_table[0]
	spamwriter = csv.writer(outf)
#	spamwriter.writerow(header)
	for row in hit_table:
		row_id = row[12] + row[22] + row[23]
		if (row_id not in spacer_mapped_dict):
			spacer_mapped_dict[row_id] = [row]
		else:
			spacer_mapped_dict[row_id].append(row)

	for mapped_hits in spacer_mapped_dict:
		i = 0
		existing_spacers = set()
		while (i < len(spacer_mapped_dict[mapped_hits])):
			row_spacer_start = spacer_mapped_dict[mapped_hits][i][8]
			row_spacer_end = spacer_mapped_dict[mapped_hits][i][9]
			index = 0
			for spacer in existing_spacers:
				spacer_start = int(float(spacer[0]))
				spacer_end = int(float(spacer[1]))
				if ((spacer_start <= int(float(row_spacer_start)) <= spacer_end or spacer_start <= int(float(row_spacer_end)) <= spacer_end or (spacer_start <= int(float(row_spacer_start)) and spacer_end >= (int(float(row_spacer_end)))) or (spacer_start >= int(float(row_spacer_start)) and spacer_end <= (int(float(row_spacer_end)))))):
					index = 1
			if (index == 0):
				existing_spacers.add((row_spacer_start, row_spacer_end))
				spamwriter.writerow(spacer_mapped_dict[mapped_hits][i])
			i += 1
	outf.close()
	return 0

# take only entries where two of more spacers map to the same sequence. Write these to a seperate file.
def two_or_more(spacer_url):
	with open (spacer_url, "r") as csvfile:
		hit_table = list(csv.reader(csvfile))

	# need to recognise when a specific array has already been detected.
	'''
	hit_dict = {}
	for hit in hit_table:	
		hit_id = hit[0].split("|")
		hit_id = hit_id[0] + hit[1] + hit[8] + hit[9] + hit[18] + hit[19] # genome_id + phage + arr_start + arr_end
		
		# remove identical targets from the same array
		if (hit_id not in hit_dict):
			hit_dict[hit_id] = hit 
	
	hit_table = list(hit_dict.values())
	'''
	hit_dict = {}
	for hit in hit_table:	
		hit_id = hit[0].split("|")
		hit_id = hit_id[0] + hit[1] + hit[18] + hit[19] # genome_id + phage + arr_start + arr_end
		
		# genome_id
		if (hit_id not in hit_dict):
			hit_dict[hit_id] = [hit] 	
		else:
			hit_dict[hit_id].append(hit)
				


	arrays = hit_dict.values()
	ret_out = open(spacer_url + "_2_or_more_hits.csv", "w")
	spam_writer = csv.writer(ret_out)

	for array in arrays:
		if (len(array) >= 2):
			for row in array:
				spam_writer.writerow(row)
	ret_out.close()
	return 0

two_or_more(sys.argv[1])




