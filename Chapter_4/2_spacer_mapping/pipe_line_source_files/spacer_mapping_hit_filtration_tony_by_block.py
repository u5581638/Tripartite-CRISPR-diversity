# script to perform spacer filtration by block prior to concatenation

import sys
import csv

spacer_only_hits_handle = open(sys.argv[1], "r")
spacer_dr_hits_handle = open(sys.argv[2], "r")
offset = int(sys.argv[4])
match_error_offset = 2
output_handle = open(sys.argv[3] + "_drs_filtered.csv", "w") # may want to make this a seperate argument

spacer_only_hits = csv.reader(spacer_only_hits_handle)
spacer_dr_hits = csv.reader(spacer_dr_hits_handle)
output_table = csv.writer(output_handle)


# strategy:
# if hit in spacer only and not in dr_spacer -> leave it to stand.
# if hit in dr and not in spacer -> remove
# if no identity difference between dr_spacer and spacer only -> remove
# if a full match to spacer yet only partial match to dr_spacer
# need to consider partial matches!!!

dr_spacer_dict = {}



for dr_spacer in spacer_dr_hits:
	dict_id = dr_spacer[0] + dr_spacer[1] # spacer_query + spacer_hit_genome. Don't include the spacer coordinates as the possibility of inexact matches mean this should be done via algebra later.
	if dict_id not in dr_spacer_dict:
		dr_spacer_dict[dict_id] = [dr_spacer]
	else:
		dr_spacer_dict[dict_id].append(dr_spacer)	

# some of this algebra is wrong!!
print("Yay!!")
for spacer in spacer_only_hits:
	spacer_dict_id = spacer[0] + spacer[1]
	if not (spacer_dict_id in dr_spacer_dict):
		output_table.writerow(spacer)

spacer_dr_hits_handle.close()
spacer_only_hits_handle.close()
output_handle.close()