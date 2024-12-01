# functions to filter spacers based on a set based approach
# No spacers should be encoded at the same position in a genome. 
# If a spacer is detected which is encoded at the same point, then assign the same coordinates as the first overlapping spacer
#
import csv 
import sys

# standardise the start/end coordinates across array predictions with different tools
def standardise(contig):
	existing_spacers = set()
	for row in contig:
		row_spacer_start = row[22]
		row_spacer_end = row[23]
		index = 0
		for spacer in existing_spacers:
			spacer_start = int(float(spacer[0]))
			spacer_end = int(float(spacer[1]))
			if ((spacer_start <= int(float(row_spacer_start)) <= spacer_end or spacer_start <= int(float(row_spacer_end)) <= spacer_end or (spacer_start <= int(float(row_spacer_start)) and spacer_end >= (int(float(row_spacer_end)))) or (row_spacer_start >= int(float(row_spacer_start)) and spacer_end <= (int(float(row_spacer_end)))))):
				row[22] = spacer_start
				row[23] = spacer_end
				index = 1
			if (index == 1):
				existing_spacers.add((row_spacer_start, row_spacer_end))	

# remove spacers with the same overlapping distance.
def same_overlap (hit_table_url, out_url):
	with open(hit_table_url, "r") as csvfile:
		hit_table = list(csv.reader(csvfile))
	# probably want to save and remove header!!
	contig_dict = {}
	outf = open(out_url,"w")
	header = hit_table[0]
	spamwriter = csv.writer(outf)
	spamwriter.writerow(header)
	for row in hit_table[1:]:
		if row[12] not in contig_dict:
			contig_dict[row[12]] = [row]
		else:
			contig_dict[row[12]].append(row)
	new_arrs = []

	for contig in contig_dict:
		i = 0
		existing_spacers = set()
		while (i < len(contig_dict[contig])):
			row_spacer_start = contig_dict[contig][i][22]
			row_spacer_end = contig_dict[contig][i][23]
			index = 0
			for spacer in existing_spacers:
				spacer_start = int(float(spacer[0]))
				spacer_end = int(float(spacer[1]))
			#	print(spacer_start)
			#	print(spacer_end)
			#	print(row_spacer_start)
			#	print(row_spacer_end)
				if ((spacer_start <= int(float(row_spacer_start)) <= spacer_end or spacer_start <= int(float(row_spacer_end)) <= spacer_end or (spacer_start <= int(float(row_spacer_start)) and spacer_end >= (int(float(row_spacer_end)))) or (spacer_start >= int(float(row_spacer_start)) and spacer_end <= (int(float(row_spacer_end)))))):
				#	print("Yay")
					contig_dict[contig][i][22] = spacer_start
					contig_dict[contig][i][23] = spacer_end
					index = 1
			if (index == 0):
				existing_spacers.add((row_spacer_start, row_spacer_end))

			i += 1
	for contig in contig_dict:
		for row in contig_dict[contig]:
			spamwriter.writerow(row)
	outf.close()
	return 0

# eliminate spacers with overlapping coordinates:
def non_redundant_hits (mapped_spacer_table_url, out_url):
	spacer_dict = {}
	with open(mapped_spacer_table_url, "r") as csvfile:
		hit_table = list(csv.reader(csvfile))
	outf = open(out_url,"w")
	spamwriter = csv.writer(outf)
	header = hit_table[0]
	spamwriter.writerow(header)
	for row in hit_table[1:]:
		genome_id = row[12]
		mapped_start = row[8]
		mapped_end = row[9]
		spacer_start = row[22]
		spacer_end = row[23]
		spacer_id = genome_id + mapped_start + mapped_end + spacer_start + spacer_end
		if (spacer_id not in spacer_dict):
			spacer_dict[spacer_id] = row 
	for my_id in spacer_dict:
		spamwriter.writerow(spacer_dict[my_id])
	outf.close()
	return 0


