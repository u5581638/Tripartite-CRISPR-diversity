# table_merger

# script to merge table of co_occurrance scores with corresponding abundance and distance scores
import csv 
import sys

# INPUT: co-occurrance table for each representative protein sequence
with open(sys.argv[1], "r") as csvfile:
	hit_table = list(csv.reader(csvfile))
# INPUT: distant/abundance table for each representative protein sequence
with open(sys.argv[2], "r") as csvfile2:
	hit_table2 = list(csv.reader(csvfile2))

hit_dict = {}

for row in hit_table2:
	if (row[0] not in hit_table2):
		hit_dict[row[0]] = row
	else:
		print("duplicate!!!!!")

ret_out = open(sys.argv[1] + "_corab_added.csv", "w")
spam_writer = csv.writer(ret_out)
for hit in hit_table:
	if (hit[0] in hit_dict):
		hit.append(hit_dict[hit[0]][1])
		my_hit = hit
		spam_writer.writerow(my_hit)
	else:
		print(hit[0])

ret_out.close()