# remove kmers matches with a dr repeat match within ~150-300bp
# Approach:
# 1. Load the mapped kmer drs.
# make non redundant. Dictionalise on the mapped sequence
# 2. Load the mapped spacers.
# iterate through checking for matches to the same genome.
# remove all instances of spacer matches where a match is found within 150-300bp on one spacer match (be strict)

import sys
import csv

spacer_open = open(sys.argv[1],"r")
spacer_kmers = list(csv.reader(spacer_open))
spacer_open.close()

dr_open = open(sys.argv[2],"r")
dr_kmers = csv.reader(dr_open)
next(dr_kmers)
print(sys.argv[1])
# first deduplicate
dr_dict = {}
for dr in dr_kmers:
	if (dr[1],dr[8],dr[9]) not in dr_dict:
		dr_dict[(dr[1],dr[8],dr[9])] = dr
dr_open.close()
dr_kmers = dr_dict.values()
# match on phage id
dr_dict = {}
for dr in dr_kmers:
	if dr[1] not in dr_dict:
		dr_dict[dr[1]] = [dr]
	else:
		dr_dict[dr[1]].append(dr)

spacer_dict = {}
for spacer in spacer_kmers:
	if (spacer[1] not in spacer_dict):
		spacer_dict[spacer[1]] = [spacer]
	else:
		spacer_dict[spacer[1]].append(spacer)

ret_spacers = []
# compare by phage_id. Exclude from return list if the target distances (spacer-drs) < 150bp. May be a high number of false negatives.
ret_out = open(sys.argv[3],"a")
spamwriter = csv.writer(ret_out)
for mapped_phage in spacer_dict:
	if (mapped_phage in dr_dict):
		dr_coords = []
		# Is [1:] here really needed?
		for dr in dr_dict[mapped_phage]:
			if(dr[8] == "Mapped_start_site"):
				continue
			dr_coords.append((int(dr[8]) + int(dr[9]))/2)

		switch = 0
		# could this be spacer_dict[mapped_phage]
		for spacer in spacer_dict[mapped_phage]:
			# extra header seems to be getting in
			if (spacer[8] == "Mapped_start_site"):
				continue
			target_coord = (int(spacer[8]) + int(spacer[9])) /2
			for dr_coord in dr_coords:
				if (abs(dr_coord - target_coord) < 150):
					switch = 1
					break
			if (switch == 1):
				break

		if (switch == 1):
			continue
		else:
			for row in spacer_dict[mapped_phage]:
				spamwriter.writerow(row)
	else:
		for row in spacer_dict[mapped_phage]:
			spamwriter.writerow(row)

ret_out.close()
