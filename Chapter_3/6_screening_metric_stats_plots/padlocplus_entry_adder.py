# pfam_entry_adder

# add DEFLOC scores to csv table containing co-occurrance scores for each representative protein sequence
from Bio import SeqIO
import sys
import csv

def row_sorter(x):
	return float(x[3])

# Input: table of co-occurrance scores
with open(sys.argv[1],"r") as csvfile2:
	co_occurrance_table = list(csv.reader(csvfile2))

# Input: Table containing DEFLOC predictions for each representative protein sequence.
with open(sys.argv[2],"r") as csvfile:
	hit_table = list(csv.reader(csvfile))


hit_dict = {}

for row in hit_table:
	if (row[0] not in hit_dict):
		hit_dict[row[0]] = [row]
	else:
		hit_dict[row[0]].append(row)
print(hit_dict.keys())
for row in hit_dict:
	hit_dict[row].sort(key=row_sorter,reverse=True)

ret_out = open(sys.argv[1] + "_padloc_added.csv","w")
spamwriter = csv.writer(ret_out)

for row in co_occurrance_table:
	if (row[0] in hit_dict):
		row.append(hit_dict[row[0]][0][1])
		row.append(hit_dict[row[0]][0][2])
		row.append(hit_dict[row[0]][0][3])
	spamwriter.writerow(row)

ret_out.close()

