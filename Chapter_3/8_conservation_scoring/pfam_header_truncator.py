import sys
import csv

# script to truncate headers in conservation matrix
# INPUT: conservation matrix
# i.e. conservation_matrix_pfam_500_all_no_CRISPR_merged.csv
# OUTPUT: conservation with shortened headers
# i.e. conservation_matrix_pfam_500_all_no_CRISPR_merged.csv_htrun.csv
# SHELL: python3 pfam_header_truncator.py conservation_matrix_pfam_500_all_no_CRISPR_merged.csv
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
