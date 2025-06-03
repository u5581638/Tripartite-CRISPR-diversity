# program to call mmseqs to recluster using a ram disk to increase the speed of writing to file
# This proccess eliminates near identical sequences from each cluster.

from Bio import SeqIO
import sys
import asyncio

# first need to retrieve clusters and put into list of lists

# takes only the first element in each cluster. Should sort first
def declusterer (reclusters):
	ret_cluster_list = []
	the_clusters = reclusters[1:]
	# cluster lacks a header (size > 1)
	cluster_head = False
	for sequ in reclusters:
		if(cluster_head == True and sequ.seq != ''):
			ret_cluster_list.append(sequ)
		if (sequ.seq == ''):
			cluster_head = True
		else:
			cluster_head = False	
	ret_cluster_list.sort(key=len, reverse=True)
	return ret_cluster_list


def mmseqs_clustering_subroutine (cluster):
	ret_Seq = []
	SeqIO.write(cluster, "cluster_file.fasta", "fasta") # must write to ramdisk dir
	subprocess.run("mmseqs", "createdb", "cluster_file.fasta", "/mnt/e/cluster_DB", "--dbtype", "1")
	subprocess.run("mkdir", "/mnt/e/tmp_cluster")
	subprocess.run("mmseqs", "cluster", "/mnt/e/cluster_DB", "/mnt/e/cluster_DB_clustered", "/mnt/e/tmp_cluster", "--min-seq-id", "0.98", "-s", "4", "cluster-mode", "2") # sensitivity shouldn't matter here hence the greedy clustering rationale. Make a note of this however.
	subprocess.run("mmseqs", "createseqfiledb", "/mnt/e/cluster_DB", "/mnt/e/cluster_DB_clustered", "/mnt/e/cluster_DB_seq")
	subprocess.run("mmseqs", "/mnt/e/cluster_DB", "/mnt/e/cluster_DB", "/mnt/e/cluster_DB_seq", "/mnt/e/cluster_DB_sequence.fasta")
	new_cluster = list(SeqIO.parse("/mnt/e/cluster_DB_sequence.fasta", "fasta"))
	# should only proccess the first sequence in list
	ret_Seq.extend(declusterer(new_cluster)) # add all entries except for the cluster head which is the first line. Sequences should already be sorted by length. Check this
	SeqIO.write(ret_Seq, "declustered_" + sys.argv[1], "fasta")
	subprocess.run("rm", "-r", "tmp_cluster")
	subprocess.run("rm", "*") # remove all files from the ramdisk
	ret_Seq = []
	return 0

# INPUT: FASTA file containing protein sequence clusters sorted in descending order and filtered by one member > 300a in length.
# i.e. one_member_gt_300_renumbered_sorted_descending_deduplicated_flattened_all_sequences.fasta_p_formatted.fastanot_in_crisprs.fasta_labelled.fasta
# OUTPUT: FASTA file containing protein sequence clusters with members in each cluster > 98% identical removed
# i.e. declustered_one_member_gt_300_renumbered_sorted_descending_deduplicated_flattened_all_sequences.fasta_p_formatted.fastanot_in_crisprs.fasta_labelled.fasta
# SHELL: python3 ramdisk_reclustering_optimised.py one_member_gt_300_renumbered_sorted_descending_deduplicated_flattened_all_sequences.fasta_p_formatted.fastanot_in_crisprs.fasta_labelled.fasta
# Note this program must be run using a ramdisk of appropiate size to support the memory requirements of many concurrent mmseqs-based clustering routines

sequences = SeqIO.parse(sys.argv[1], "fasta")

sequence_clusters = []
cluster = []
past_cluster_name = "random pre_set"
out_file = open("declustered_" + sys.argv[1], "a")

for sequ in sequences:
	this_sequ = sequ.description.split("|cluster_")
	this_sequ = this_sequ[-1]
	print(this_sequ)
	if (this_sequ != past_cluster_name):
		mmseqs_clustering_subroutine(cluster) # will initially go through and append first cluster
		past_cluster_name = this_sequ
	cluster.append(sequ)
sequence_clusters.append(cluster)
print("ready")
# then need to: 
# 1. Write each cluster to a RAMdisk
# 2. Via subproccess.run():
# 3. run createdb 
# 4. running clustering
# 5. flatten file
# 6. rm all intermediate files
# 7. load cluster back into python and remove the top cluster header before writing to a return file.


print("Clustering complete!")
	



