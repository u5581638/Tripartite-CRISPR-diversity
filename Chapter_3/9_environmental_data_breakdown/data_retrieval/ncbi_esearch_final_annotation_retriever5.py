# ncbi_esearch_final_annotation_retriever
# script to search for and retrieve species names by ncbi identifier using entrez

import sys
import csv
import re
from Bio import Entrez
# need to:
# 1. use re to extract ncbi only identifers post ncbi_table based extraction run
# 2. Enable flexible resizeing of the input string.
# output a table giving at least the ncbi species names.
def row_sorter(x):
	return x[0]

# SHELL: see large_scale_data_block_species_rncbi_retriever6_annotation_merged.sh
# This script was run recursively using half the batch_number with each iteration.
# INPUT: 1. table containing NCBI identifiers compatible with entrez
#		 2. number of identifers to query at once (int)
# OUTPUT: 1. table containing a list of matching ncbi identifiers
#		  2. List of unmatched identifiers
#		  1. list of incompatible IDs
#		   

batch_num = int(sys.argv[2])

# table/list of ncbi identifiers in csv format
with open (sys.argv[1],"r") as csvfile:
	hit_table = csv.reader(csvfile)
	hit_dict = {}
	for hit in hit_table:
		hit_dict[hit[0]] = hit
	hit_table = list(hit_dict.values())

hit_table.sort(key=row_sorter)
ret_list = []
ele_str = ""
# OUTPUT:
# 1. retrieved records
ret_out = open(sys.argv[1] + "_ncbi.csv","a")
spamwriter = csv.writer(ret_out)
# remaining identifers where retrieval failed
remaining_out = open(sys.argv[1] + "_rncbi.csv","a")
remwriter = csv.writer(remaining_out)
# Other identifiers where retrieval is not possible.
other_out = open(sys.argv[1] + "_oncbi.csv","a")
otherwriter = csv.writer(other_out)
max_ele_size = 0
id_block_dict = {}
for  identifier in hit_table:
#	print(re.search('.*\.[0-9]',identifier[0]))
#	print(not re.search('Contig.*$',identifier[0]))
#	print(not re.search('NODE.*$',identifier[0]))
	block_id = identifier[0].split("|")[0]

	identifier = identifier[0].split("|")[1]
	id_block_dict[identifier] = block_id

	if (re.search('.*\.[0-9]',identifier) and not re.search('Contig.*$',identifier) and not re.search('NODE.*$',identifier) and not re.search('_',identifier)):
		ele_str = ele_str + identifier + ","
		max_ele_size += 1
		if (max_ele_size > batch_num ):
			ret_list.append(ele_str[:-1])
			max_ele_size = 0
			ele_str = ""
	else:
		otherwriter.writerow([block_id + "|" + identifier])
other_out.close()
ret_list.append(ele_str)

for query in ret_list:
	if (query == ""):
		continue

	Entrez.email = "u5581638@anu.edu.au"
	handle = Entrez.esearch(db="nucleotide", term=query,retmax=str(batch_num + 10))
	try:
		records = Entrez.read(handle)
	except:
		rem_queries = query.split(",")
	#	print(rem_queries)
		for q in rem_queries:
			remwriter.writerow([id_block_dict[q] + "|" + q])
		continue	
	gid_str = ""
	for rec in records["IdList"]:
		gid_str = gid_str + rec + ","
	gid_str = gid_str[:-2]

	
	Entrez.email = "u5581638@anu.edu.au"
	handle = Entrez.esummary(db="nucleotide",id=gid_str)
	try:
		records = Entrez.read(handle)
		for record in records:
			spamwriter.writerow([record['AccessionVersion'],record['Title']])
	except:
		rem_queries = query.split(",")
		for q in rem_queries:
			remwriter.writerow([id_block_dict[q] + "|" + q])
		print(query)
		print("no_compatible_records!")	
		

ret_out.close()
remaining_out.close()





