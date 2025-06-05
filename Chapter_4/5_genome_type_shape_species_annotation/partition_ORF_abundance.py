import csv
import sys
import os

# INPUT: 1. table showing partition numbers for each node from leiden clustering
#		 i.e. cas12a_host_only_partition.csv_top10.csv
#		 2. protein annotation table for each CRISPR-Cas subtype or the corresponding mapped target sequences
#		 i.e. cas12a.fasta_all_hits.csv_genomes.fasta_annotations.csv
# OUTPUT: set of table files (in .csv format) grouping the annotations by partitioned cluster. Each cluster is a seperate file. These are in a "results/" folder
# SHELL: python3 partition_ORF_abundance.py cas12a_host_only_partition.csv_top10.csv cas12a.fasta_all_hits.csv_genomes.fasta_annotations.csv 

# open the partition table
with open(sys.argv[1]) as csvfile:
	partition_table = list(csv.reader(csvfile))

# Open annotation table. For phages, I'd need to use mapped phage ORFs. This may not be feasible given the sheer number of phages. Need to do some kind of parametric screen. (i.e. minimum cluster size > 20)
with open(sys.argv[2]) as csvfile:
	annotation_table = list(csv.reader(csvfile))

# need to group each partition then retrieve the corresponding annotations. Each group of ORFs should be written to file then counted by a seperate script.

partition_dict = {}

for row in partition_table[1:]:
	if row[2] not in partition_dict:
		partition_dict[row[2]] = [row[1]]
	else:
		partition_dict[row[2]].append(row[1])

annotation_dict = {}

for row in annotation_table[1:]:
	if (row[0] not in annotation_dict):
		annotation_dict[row[0]] = [row]
	else:
		annotation_dict[row[0]].append(row)

os.mkdir("results/")
for group in partition_dict:
	partition = partition_dict[group]
	ret_out = open("results/" + "cluster_" + group + ".csv", "w")
	spamwriter = csv.writer(ret_out)
	for row in partition:
		if (row in annotation_dict):
			for r in annotation_dict[row]:
				spamwriter.writerow(r)
		else:
			print(row)
	ret_out.close()



