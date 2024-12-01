# program to merge the color and partition tables (w/ color column)
import sys
import csv


# table containing list of sequences grouped into clusters by leiben partitioning
with open(sys.argv[1],"r") as csvfile:
	partition_table = list(csv.reader(csvfile))

# table of cluster specific colours
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

