import csv
import sys

# filter genomes with homology to RT sequence (by BLASTp/BLASTx) with those without.
with open(sys.argv[1],"r") as csvfile:
	blast_table = list(csv.reader(csvfile))

blast_dict = {}
for row in blast_table:
	blast_dict[row[1]] = row 
blast_table = list(blast_dict.values())
with open(sys.argv[2],"r") as csvfile:
	mapping_table = list(csv.reader(csvfile))

mapped_dict = {}

for row in mapping_table[1:]:
	my_id = row[0].split("|") [0]
	if my_id not in mapped_dict:
		mapped_dict[my_id] = [row]
	else:
		mapped_dict[my_id].append(row)


rt_row_dict = {}
ret_out = open(sys.argv[3],"w")
spamwriter = csv.writer(ret_out)
spamwriter.writerow(mapping_table[0])
for row in blast_table:
	if row[1] in mapped_dict:
		matches = mapped_dict[row[1]]
		rt_row_dict[row[1]] = mapped_dict[row[1]]
		for r in matches:
			spamwriter.writerow(r)

ret_out.close()
ret_out4 = open(sys.argv[4],"w")
spamwriter2 = csv.writer(ret_out4)
spamwriter2.writerow(mapping_table[0])
for row in mapping_table[1:]:
	if (row[0].split("|")[0] in rt_row_dict):
		spamwriter2.writerow(row)
ret_out4.close()