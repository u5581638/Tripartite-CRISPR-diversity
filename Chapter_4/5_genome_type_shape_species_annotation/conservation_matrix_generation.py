#conservation_matrix_generation
import sys
import csv

# generate a matrix of conservation scores from an input list of csv tables containing the conservation score for each gene annotation.
# INPUT: list of conservation tables for each partitioned cluster
# i.e. cluster_1.csv_anno.csv cluster_2.csv_anno.csv ... cluster_n.csv_anno.csv
# OUTPUT: matrix showing conservation scores for the same conserved protein annotations from each partitioned cluster.
# i.e. conservation_matrix_pfam_500_all.csv
# SHELL: python3 conservation_matrix_generation.py  .. etc <output_matrix_name>

i = 1
conservation_dict = {}
while (i < len(sys.argv) - 1):
	with open(sys.argv[i],"r") as csvfile:
		con_table = list(csv.reader(csvfile))
	#	print(con_table)
		# use if clause depending on whether table has a header
		if (con_table[0][0] == "Protein name"):
			conservation_dict[sys.argv[i]] = con_table[1:]
		else:
			conservation_dict[sys.argv[i]] = con_table	
	i += 1

headers = set()
# add all heading names
for my_table in conservation_dict:
	table = conservation_dict[my_table]
	for row in table:
		if(row[0] not in headers):
			headers.add(row[0])

all_headers = list(headers)
all_headers.sort()
all_headers = [""] + all_headers
ret_list = [all_headers]
# how to build each row. Header/x axis is new proteins. Y is subtype.
# need to realign each value to the header
for subtype in conservation_dict:
	retrow = [subtype.split(".")[0]] + ["0"] * (len(all_headers) - 1)
	for row in conservation_dict[subtype]:
		index = all_headers.index(row[0])
	#	print(row,subtype)
		print(index,row)
		retrow[index] = row[1]
	ret_list.append(retrow)

ret_out = open(sys.argv[len(sys.argv) - 1],"w")
spamwriter = csv.writer(ret_out)
for ret in ret_list:
	spamwriter.writerow(ret)
ret_out.close()