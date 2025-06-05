import sys
import csv
import re

# script to filter unassociated/singleton sequences from the interaction table.
# INPUT: host-phage interaction table (from homolog_node_constructor_numbering_tabulation.py)
#		i.e cas13b_host_phage_interaction_table.csv
# OUTPUT: host-phage interaction table with singleton entries removed.
#		i.e. cas13b_host_phage_interaction_table.csv_names_kept2.csv
csvfile = open(sys.argv[1])
hp_table = csv.reader(csvfile)

ret_out = open(sys.argv[1] + "_names_kept2.csv","w")
spamwriter = csv.writer(ret_out)
spamwriter.writerow(next(hp_table))
hp_set = set()
for row in hp_table:
	host_name = row[2]
	phage_name = row[3]
	my_matchh = re.match('h[0-9]+',host_name)
	my_matchp = re.match('p[0-9]+',phage_name)
	if (my_matchh and my_matchp and ((host_name, phage_name) not in hp_set)):
		spamwriter.writerow(row)
		hp_set.add((host_name, phage_name))
ret_out.close()
csvfile.close()
