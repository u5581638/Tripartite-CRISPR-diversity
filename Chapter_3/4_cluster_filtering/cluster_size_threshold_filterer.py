# script to filter clusters of size >= 3
from Bio import SeqIO
import sys

# Input file containing clusters in FASTA format.
sequences = SeqIO.parse(sys.argv[1], "fasta")
cluster_size_cutoff = 3
cluster_size = 0
sequ = next(sequences)
this_sequ = sequ.description.split("cluster_")
this_sequ = this_sequ[-1]
previous_sequ = this_sequ
cluster_size += 1
cluster = [sequ]

out_file = open("grt_" + str(cluster_size_cutoff) + "_" + sys.argv[1], "a")
for sequ in sequences:
	this_sequ = sequ.description.split("cluster_")
	this_sequ = this_sequ[-1]
	if (this_sequ != previous_sequ):
		print(cluster_size)
		if (cluster_size >= cluster_size_cutoff):
			print("cluster size is:" + str(len(cluster)))
			SeqIO.write(cluster, out_file, "fasta")
		cluster = []
		cluster_size = 0
	cluster_size += 1
	cluster.append(sequ)
	previous_sequ = this_sequ
#print("complete!")
out_file.close()

