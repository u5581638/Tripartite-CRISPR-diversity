import sys
import csv

# script to truncate headers in conservation matrix

with open(sys.argv[1],"r") as csvfile:
	the_matrix = list(csv.reader(csvfile))

header = the_matrix[0]
new_header = []
for my_id in header:
	my_id = my_id.split(" ")
	if (len(my_id) > 1):
		new_id = my_id[0] + " " + my_id[1]
	else:
		new_id = my_id[0]
	new_header.append(new_id)

ret_out = open(sys.argv[1] + "_htrun.csv","w")
spamwriter = csv.writer(ret_out)
spamwriter.writerow(new_header)
for row in the_matrix[1:]:
	spamwriter.writerow(row)
ret_out.close()
