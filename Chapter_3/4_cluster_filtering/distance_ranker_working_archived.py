# distance_ranker
# return a ranked list of clusters and representative sequences based on the distance from the CRISPR array

from Bio import SeqIO
import sys
import re

# input set of representative protein sequences with average distance labelled in identifiers (in FASTA format)
sequences = SeqIO.parse(sys.argv[1], "fasta")

distance_pattern = re.compile("Average_dist=[0-9]*")
maximum_distance = 2000 # 2kb from CRISPR array
ret_seq = []

for sequ in sequences:
#	print(sequ.description)
	avg_distance_regex = re.search(distance_pattern, sequ.description)
	avg_distance_str = avg_distance_regex.group(0)
	avg_distance = int(avg_distance_str[13:])
	if (avg_distance < maximum_distance):
		ret_seq.append(sequ)
# returns the clusters with a distance less than the maximum distance.		
SeqIO.write(ret_seq, "ranked_" + str(maximum_distance) + sys.argv[1], "fasta")

# now just need to select the longest representative sequence from each cluster and we're done!!

def longest_sequence(sequence_list):
	sequence_length = 0
	long_sequ = sequence_list[0]
	for sequ in sequence_list:
		if (len(sequ.seq) > sequence_length):
			sequence_length = len(sequ)
			long_sequ = sequ
	return long_sequ		

representative_sequences = []
cluster_pattern = re.compile("cluster_[0-9]* CRISPR_Position")
sequ = ret_seq[0]
first_sequence_regex = re.search(cluster_pattern, sequ.description)
first_sequence_pattern = first_sequence_regex.group(0)

previous_cluster_no = first_sequence_pattern.split("cluster_")
previous_cluster_no = previous_cluster_no[1]
previous_cluster_no = previous_cluster_no.split(" CRISPR_Position")
previous_cluster_no = int(previous_cluster_no[0])

cluster = [sequ]

i = 1
while (i < len(ret_seq)):
	current_sequence_regex = re.search(cluster_pattern, ret_seq[i].description)
	current_cluster_pattern = current_sequence_regex.group(0)
	current_cluster = current_cluster_pattern.split("cluster_")
	current_cluster = current_cluster[1]
	current_cluster = current_cluster.split(" CRISPR_Position")
	current_cluster = int(current_cluster[0])
	if (current_cluster != previous_cluster_no):
		longest_seq = longest_sequence(cluster)
		representative_sequences.append(longest_seq)
		cluster = []
		cluster.append(ret_seq[i])
	else:
		cluster.append(ret_seq[i])
	i += 1
	previous_cluster_no = current_cluster
# returns the longest representative sequence within each cluster with average distance from the CRISPR array less than the maximum distance	
SeqIO.write(representative_sequences, "representative_distance_" + sys.argv[1], "fasta")			