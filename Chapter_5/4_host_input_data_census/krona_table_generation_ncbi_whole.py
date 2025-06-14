# krona_table_generation

# need to sort and create tables compatible with the krona plot

import sys
import csv

# INPUT: NCBI table containing species metadata for each subtype in csv format.
# Note: This was run for JGI/NCBI assembles using two variants of this script.
# OUTPUT: table of contig annotations usable as input to Krona.html

# i.e. ncbi_microbiome_window_only2_dedup.csv
with open(sys.argv[1],"r") as csvfile:
	raw_table = list(csv.reader(csvfile))

ncbi_dict = {}

# need to use tuple (datatype,genus,species) as tuple
# need to count number of instances and store as value
# unusual inversion of standard dict use.

for row in raw_table:
	my_key = []
	 # change this index for jgi
	if (len(row) > 2):
		gen_species = row[8].split(" ")	
		datatype = row[25]
	elif(len(row) == 2):
		gen_species = row[1].split(" ")	
		datatype = "Unclassified"
	else:
		continue
	my_key.append("NCBI")
	my_key.append(datatype)
	if (len(gen_species) > 1 and gen_species[1] != "metagenome"):
		genus = gen_species[0]
		species = gen_species[1]
		my_key.append(genus)
		my_key.append(species)
	else:
		genus = gen_species[0]
		species = ""
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

