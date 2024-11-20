import csv
import sys

with open(sys.argv[1],"r") as csvfile:
	interaction_table = list(csv.reader(csvfile))

with open(sys.argv[2],"r") as csvfile:
	host_node_table = list(csv.reader(csvfile))

# need phage_w_color.csv table from the bioinfo computer to make this work!!
with open(sys.argv[3],"r") as csvfile:
	phage_node_table = list(csv.reader(csvfile))

host_dict = {}
phage_dict = {}
for row in host_node_table:
	if (row[1] not in host_dict):
		host_dict[row[1]] = row

for row in phage_node_table:
	if (row[1] not in phage_dict):
		phage_dict[row[1]] = row

ret_out = open(sys.argv[1] + "_hp_colors_table_top10_only.csv","w")
spamwriter = csv.writer(ret_out)
spamwriter.writerow(interaction_table[0] + ["host_cluster", "host_color","phage_cluster","phage_color"])
for row in interaction_table[1:]:
	myrow = row
	if (row[0] in host_dict):
		myrow = myrow + host_dict[row[0]] [2:4]
	else:	
		myrow = myrow + ["NA","#000000"]
	#	print("not in row!!")
	if (row[4] in phage_dict):
		myrow = myrow + phage_dict[row[4]] [2:4]
	else:	
		myrow = myrow + ["NA","#000000"]
	spamwriter.writerow(myrow)
