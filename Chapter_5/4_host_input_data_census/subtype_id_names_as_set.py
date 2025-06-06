# subtype_id_names_as_set

# script to take the mapped spacer tables and extract just the genome ids for proportion calculation

# must add top layer to krona layers to be either JGI or NCBI -> this will enable merging of the two seperate plots.

import sys
import csv
from Bio import SeqIO

reader_url = open(sys.argv[1],"r")
mapping_handle = csv.reader(reader_url)

ret_set = {}
next(mapping_handle)
for row in mapping_handle:
	my_id = row[0].split("|")[0]
	if my_id not in ret_set:
		ret_set[my_id] = my_id

reader_url.close()
ret_out = open(sys.argv[2],"w")
spamwriter = csv.writer(ret_out)

for row in ret_set.keys():
	spamwriter.writerow([row])

ret_out.close()
