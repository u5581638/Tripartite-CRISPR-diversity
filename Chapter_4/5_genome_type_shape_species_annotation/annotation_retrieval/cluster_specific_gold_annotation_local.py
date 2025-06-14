# retrieve the subset of gold annotations specific to a subset of host genomes

# Need to potentially index the csv files to enable fast lookup.

# THis will be different for host vs phage genomes because local genomes do not have a block identifier!!


# 2nd best strategy is to sort entries by block so the minimal number of load operations can be used
# 1. Load interaction table
# 2. sort copy of interaction table
# 3. Remove redundant entries
# 4. Lookup relevant labelled genomes file

import csv
import sys

# INPUT: 1. Spacer-mapping interaction table specific to each CRISPR-Cas subtype
#		 2. Gold annotation table for CRISPR 40kb DNA windows
# OUTPUT: 1. Table of annotations matched to just the host-encoded CRISPR-Cas subtypes.

with open(sys.argv[1],"r") as csvfile:
	interaction_table = list(csv.reader(csvfile))

csvfile2 = open(sys.argv[2], "r")
gold_table = csv.reader(csvfile2)

gold_dict = {}

for row in gold_table:
	if (row[0].split("_")[0] not in gold_dict):
		gold_dict[row[0].split("_")[0]] = row

csvfile2.close()
my_out = open(sys.argv[3],"a")
spamwriter = csv.writer(my_out)

known_rows = set()
for row in interaction_table:
	my_row = row[0].split("_") [0]
	if (my_row in gold_dict and my_row not in known_rows):
		known_rows.add(my_row)
		spamwriter.writerow(gold_dict[my_row])
	else:
		print("missing_identifer!!")
		print(row)
my_out.close()


