# krona_table_generation

# need to sort and create tables compatible with the krona plot

import sys
import csv

# INPUT: JGI table containing species metadata in csv format.
with open(sys.argv[1],"r") as csvfile:
	raw_table = list(csv.reader(csvfile))

ncbi_dict = {}


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

