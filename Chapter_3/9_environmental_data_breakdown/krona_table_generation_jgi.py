# krona_table_generation

# need to sort and create tables compatible with the krona plot

import sys
import csv

# INPUT: JGI table containing species metadata in csv format.
# Note: This was run for JGI/NCBI assembles using two variants of this script.
# i.e. jgI_microbiome_window_only2.csv

with open(sys.argv[1],"r") as csvfile:
	raw_table = list(csv.reader(csvfile))

ncbi_dict = {}

# OUTPUT: table of contig annotations usuable as input to Krona.html
# SHELL: krona_table_generation_jgi.py jgI_microbiome_window_only2.csv <output_file_name>
for row in raw_table:
	my_key = []
#	datatype = row[24] # change this index for jgi
#	gen_species = row[7].split(" ")	
#	my_key.append(datatype)
	my_class = row[7]
	order = row[8]
	genus = row[9]
	species = row[10]
	my_key.append(my_class)
	my_key.append(order)
	my_key.append(genus)
	my_key.append(species)
	my_key = tuple(my_key)
	if (my_key not in ncbi_dict):
			ncbi_dict[my_key] = 1
	else:
		ncbi_dict[my_key] += 1

ret_out = open(sys.argv[2],"w")
spamwriter = csv.writer(ret_out)
for tkey in ncbi_dict:
	spamwriter.writerow(list(tkey) + [str(ncbi_dict[tkey])])
ret_out.close()				

