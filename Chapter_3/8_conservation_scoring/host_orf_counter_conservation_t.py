# script to count the number of orf annotations then to compute the conservation scores

import sys
import csv

# INPUT: 1 protein annotation table generated for each CRISPR-Cas subtype (by annotation_parallelisation_function_SQL11.py - see chapter 4)
# i.e. cas13b.fasta_all_hits.csv_genomes.fasta_annotations_all.csv
# OUTPUT: list of each protein annotation with it's conservation score annotated.
# i.e. cas13b.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_conservation.csv
# SHELL: host_orf_counter_conservation_t.py cas13b.fasta_all_hits.csv_genomes.fasta_annotations_all.csv
with open(sys.argv[1],"r") as csvfile:
	proteins = list(csv.reader(csvfile))[1:]

prot_dict = {}

for row in proteins:
	prot_dict[row[1]] = row

proteins = list(prot_dict.values())
prot_index = {}
index = 19
global total_genomes
total_prots = []
for row in proteins:
	total_prots.append(row[0])
	my_row = row[index].split(" ")
	my_row = my_row[0]
	if (my_row not in prot_index):
		prot_index[my_row] = 1
	else:
		prot_index[my_row] += 1	 

total_prots = list(set(total_prots))
total_genomes = len(total_prots)
def sorter(x):
	return int(x[1])

def conservation_proportion(x):
	return (x[0], int(x[1]) / total_genomes)

prot_index = list(prot_index.items())
prot_index.sort(key=sorter,reverse=True)
prot_index = map(conservation_proportion,prot_index)
my_out = open(sys.argv[2],"w")
spamwriter = csv.writer(my_out)
max_rows = 500
i =0
headers = ["Protein name"]
numbers = ["Conservation %"]
spamwriter.writerow(["Protein name","Conservation"])
for rw in prot_index:
	if (i == max_rows):
		break
	spamwriter.writerow(list(rw))
	i += 1
	headers.append(rw[0])
	numbers.append(rw[1])

#spamwriter.writerow(headers)
#spamwriter.writerow(numbers)	
my_out.close()



