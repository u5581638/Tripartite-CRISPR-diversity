# program to relabel clusters after declustering

from Bio import SeqIO
import sys

# files containing labelled clusters in FASTA format.
sequences = SeqIO.parse(sys.argv[1], "fasta")
cluster_number = 1
sequ = next(sequences)
sequ.description = sequ.description + "|cluster_" + str(cluster_number)
this_sequ = sequ.description.split("cluster_")
this_sequ = this_sequ[-1]
previous_sequ = this_sequ
cluster = [sequ]

outfile = open("relabelled_" + sys.argv[1], "a")
for sequ in sequences:
	this_sequ = sequ.description.split("cluster_")
	this_sequ = this_sequ[-1]
	if (this_sequ != previous_sequ):
		SeqIO.write(cluster, outfile , "fasta")
		cluster = []
		cluster_number += 1
	sequ.description = sequ.description + "|cluster_" + str(cluster_number)	
	cluster.append(sequ)
	previous_sequ = this_sequ
SeqIO.write(cluster, outfile, "fasta")	
outfile.close()
print("total_number of clusters is: " + str(cluster_number))
