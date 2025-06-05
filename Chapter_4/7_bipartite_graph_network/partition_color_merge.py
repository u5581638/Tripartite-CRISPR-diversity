# program to merge the color and partition tables (w/ color column)
import sys
import csv

# INPUT: 1. table containing list of sequences grouped into clusters by leiben partitioning (from "retrieve_graph_partitions.r")
#		 i.e. cas12a_host_only_partition_top10.csv
#		 2. table of partition specific colours (from "vCONTACT_only_graph_generation2_partition_NA_top10.r")
#		 i.e cas12a_host_only_partition_top10.csv_color_table.csv
# OUTPUT: table assigning a colour and partition to all nodes
#		 i.e. cas12a_host_only_partition_top10_color_only.csv
# SHELL: python3 partition_color_merge.py cas12a_host_only_partition_top10.csv cas12a_host_only_partition_top10.csv_color_table.csv cas12a_host_only_partition_top10_color_only.csv
# table containing list of sequences grouped into clusters by leiden partitioning
with open(sys.argv[1],"r") as csvfile:
	partition_table = list(csv.reader(csvfile))

# table of partition specific colours
with open(sys.argv[2],"r") as csvfile:
	color_table = list(csv.reader(csvfile))

color_dict = {}

for row in color_table[1:]:
	if row[1] not in color_dict:
		color_dict[row[1]] = row[4]
	else:
		color_dict[row[1]] = row[4]

ret_out = open(sys.argv[3],"w")
spamwriter = csv.writer(ret_out)
spamwriter.writerow(["","node","component","color"])
for row in partition_table[1:]:
	myrow=row
	if row[2] in color_dict:
		myrow.append(color_dict[row[2]])
		spamwriter.writerow(myrow)
	else:
		print("error!!!!!")	

ret_out.close()

