# distance calculator
# Find the minimum distance between the array and the starting point of the arrays
# this requires computing the abolsute positions of the proteins from the CRISPR array
# the abolsute position will be the position given by prodigal/genemark + starting array position
# if strand is negative the starting position should have a higer score than the negative position
# do this for each protein in each cluster, then average the results for each cluster - this is now the primary purpose of the program is the previous program finds the minimal distance and  starting point of the arrays.

from Bio import SeqIO
import sys
import re
from statistics import mean
from statistics import stdev

# file containing filtered clusters with distances labelled in the identifiers in FASTA format.
sequences = list(SeqIO.parse(sys.argv[1], "fasta"))
i = 0
cluster_pattern = re.compile("cluster_[0-9]* CRISPR_Position")
first_sequence_regex = re.search(cluster_pattern, sequences[i].description)
first_sequence_pattern = first_sequence_regex.group(0)
print(first_sequence_pattern)
previous_cluster_no = first_sequence_pattern.split("cluster_")
previous_cluster_no = previous_cluster_no[1]
previous_cluster_no = previous_cluster_no.split(" CRISPR_Position")
previous_cluster_no = int(previous_cluster_no[0])
 # "extracted cluster number!"
distances = []

distance_pattern = re.compile("Distance=[0-9]*")
first_sequence_regex_distance = re.search(distance_pattern, sequences[i].description)
first_sequence_distance = first_sequence_regex_distance.group(0)
first_distance = int(first_sequence_distance[9:])

distances.append(first_distance)
sequence_indexes = []
sequence_indexes.append(i) 
i = 1
while (i < len(sequences)):
	current_cluster_sequence = re.search(cluster_pattern, sequences[i].description)
	current_cluster_pattern = current_cluster_sequence.group(0)
#	print(current_cluster_pattern)
	current_cluster = current_cluster_pattern.split("cluster_")
	current_cluster = current_cluster[1]
	current_cluster = current_cluster.split(" CRISPR_Position")
	current_cluster = int(current_cluster[0])
	if (current_cluster != previous_cluster_no): # need to execute the averaging subfunction and update all sequences concerned with this subdistance.
		avg_distance = mean(distances)
		stdev_score = stdev(distances)
		print(distances)
		k = 0
		while (k < len(sequence_indexes)):
			sequences[sequence_indexes[k]].description += " " + "Average_dist=" + str(avg_distance) + " " + "Stdev=" + str(stdev_score)
			k += 1
		sequence_indexes = []
		distances = []
		seq_regex_distance = re.search(distance_pattern, sequences[i].description)
		seq_distance = seq_regex_distance.group(0)
		distance = int(seq_distance[9:])
		distances.append(distance)
		sequence_indexes.append(i)	
	else:
		seq_regex_distance = re.search(distance_pattern, sequences[i].description)
		seq_distance = seq_regex_distance.group(0)
		distance = int(seq_distance[9:])
		distances.append(distance)
		sequence_indexes.append(i)
	previous_cluster_no = current_cluster	
	i += 1

avg_distance = mean(distances)
stdev_score = stdev(distances)
print(distances)
k = 0
while (k < len(sequence_indexes)):
	sequences[sequence_indexes[k]].description += " " + "Average_dist=" + str(avg_distance) + " " + "Stdev=" + str(stdev_score)
	k += 1	

SeqIO.write(sequences, "avg_dist_" + sys.argv[1], "fasta") # output sequences with all distances appended


