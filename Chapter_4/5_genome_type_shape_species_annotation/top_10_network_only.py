import csv
import sys

# INPUT: table containing partitions delineated by leiden clustering.
# 		 i.e. cas12a_host_only_partition.csv
# OUTPUT: table containing only the top 10 largest partitions
#		 i.e. cas12a_host_only_partition_top10.csv
# SHELL: python3 top_10_network_only.py <input_file.csv> <output_file.csv>


with open(sys.argv[1],"r") as csvfile:
	hit_table = list(csv.reader(csvfile))

hit_dict = {}

for hit in hit_table[1:]:
	if (hit[2] not in hit_dict):
		hit_dict[hit[2]] = 1 
	else:
		hit_dict[hit[2]] += 1

dict_items = list(hit_dict.items())
dict_items = sorted(dict_items,key=lambda item: item[1],reverse=True)
clade_number = 0
cluster_numbers = set()
for group in dict_items:
	if (clade_number >= 10):
		break
	if (group[0] not in cluster_numbers):
		cluster_numbers.add(group[0])
	clade_number += 1

ret_out = open(sys.argv[2],"w")
spamwriter = csv.writer(ret_out)
# write header
spamwriter.writerow(hit_table[0])
for hit in hit_table[1:]:
	if (hit[2] in cluster_numbers):
		spamwriter.writerow(hit)
	else:
		my_hit=hit 
		my_hit[2] = "NA"
		spamwriter.writerow(my_hit)
ret_out.close()