# cluster sorter
# script to sort individual clusters by size, then sort clusters in descending order of total cluster size

from Bio import SeqIO
import sys
import itertools

sequences = SeqIO.parse(sys.argv[1], "fasta")

clusters = []
cluster = []
cluster_number = 0
for sequ in sequences:
	this_sequ = sequ
	this_sequ.description = this_sequ.description + "|" + "cluster_" + str(cluster_number)
	if (this_sequ.seq == ''):
		clusters.append(cluster)
		cluster = []
		cluster_number += 1
	else:
		cluster.append(this_sequ)
clusters.append(cluster)		
clusters = clusters[1:] # remove empty list element at the start of group

if (clusters[-1] == [] ):
	clusters = clusters[:-2]

print("grouping complete!")

# next sort each cluster by len of sequences.

for clust in clusters:
	clust.sort(key=len, reverse=True)
clusters.sort(key=len, reverse=True)

flattened_list = []

for clust in clusters:
	flattened_list.extend(clust)
SeqIO.write(flattened_list, "sorted_descending_" + sys.argv[1], "fasta")	