# script to filter clusters of size >= 3
from Bio import SeqIO
import sys

# INPUT: Output from cluster_length_filter after filtering members > 98% similar. A list of protein sequences labelled by cluster in FASTA format, with clusters lacking at least one sequence > 300aa removed.
# i.e. cluster_declustered_one_member_gt_300_renumbered_sorted_descending_deduplicated_flattened_all_sequences.fasta_p_formatted.fastanot_in_crisprs.fasta_labelled.fasta
# OUTPUT: FASTA file containing protein sequences labelled by cluster, with each cluster >= 3 members in size.
# i.e. grt_3_cluster_declustered_one_member_gt_300_renumbered_sorted_descending_deduplicated_flattened_all_sequences.fasta_p_formatted.fastanot_in_crisprs.fasta_labelled.fasta
# SHELL: python3 cluster_size_threshold_filterer.py cluster_declustered_one_member_gt_300_renumbered_sorted_descending_deduplicated_flattened_all_sequences.fasta_p_formatted.fastanot_in_crisprs.fasta_labelled.fasta

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

