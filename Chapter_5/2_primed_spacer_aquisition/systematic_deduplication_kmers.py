# Deduplicate mapped hits from identical or near identical coordinates.
'''
Approach:
# Two possible ways to recognise these hits:
# 1. By DR-spacer - spacer comparison. This is problematic for kmers.
# What if I do an equivalent kmer search using the direct repeats as input, then exclude based on where a direct repeat hit occurs within ~100-150bp of a spacer kmer.
# This should work, hypothetically at least.
# Will need two scripts -> one to perform the equivalent kmer search using array DRs -> This should in theory be much easier because the DR is the same.
# The second script will need to perform the exclusion
# may be work run this script on other data to ensure the removal of self-targeting CRISPR-arrays although this will likely result on the loss of some
important phenomena as well.
# 2. By Proximity to the projected spacer. This is also problematic because the coordinates of the kmer in the spacer are unknown and not regular.

Approach 2:

1. Deduplicate by distance -> for a given phage-host pair. None of the respective distances when looking at the pair should be within a certain minimum distance of each other. This includes distances calculated not just from the PPS but from other spacers as well.

2. Need to deduplicate the case of two arrays mapping to the same spacer i -> if the absolute distances between the spacers are the same then keep just one.

In both cases, either extending the length of the kmers or recording the original size of the spacers, as well as the kmerised section, will be essential to calculating where duplication is occuring. The only alternative is filtering distances over a much greater range.
'''

import csv
import sys

# need to consider the length and positions of the kmers to compute an allowed distance between them.
# Approach:
# iterate through each group.
# create a ret_dict for each group
# if the distances between the spacers exceed the allowed group size then add. This includes non-PPS distances.
# need to prioritise keeping pairs of hits form the same array/same phage.
def distance_based_deduplication_array_names(host_spacer_pairs,allowed_distance_between_pairs=50):
	ret_host_spacer_pairs = []
	# iterate through each possible host target pairing
	distance_matrix = []
	a = 0
	# dictionary to contain array identifiers.
	# iterate through each set of arrays/mapped phages:
	while (a < len(host_spacer_pairs)):
		i = 0 
		array_names = {} 
		while (i < len(host_spacer_pairs[a])):
			allowed_overlap = allowed_distance_between_pairs
			k = i + 1
			while (k < len(host_spacer_pairs[a])):			
				host_spacer_distance = abs(float(host_spacer_pairs[a][i][-2]) - float(host_spacer_pairs[a][k][-2])) # index should be the mapped distance.
				query_start = int(host_spacer_pairs[a][k][6])
				query_end = int(host_spacer_pairs[a][k][7])
				query_size = abs((query_start - query_end))
				query_length = int(host_spacer_pairs[a][k][5])
				# could change this to 20.
			#	allowed_overlap = query_length - 2 * query_size 
			#	allowed_overlap = 20
				if (host_spacer_distance < allowed_overlap): # If this happens rather than eliminating, should instead group the elements together.
				# May need to pop to prevent reconsideration of the same entries
					if ((host_spacer_pairs[a][i][0], host_spacer_pairs[a][i][-4]) in array_names):
						# add name to standardised genome_id mapp
						array_names[(host_spacer_pairs[a][k][0],host_spacer_pairs[a][k][-4])] = (host_spacer_pairs[a][i][0], host_spacer_pairs[a][i][-4])
						# standardise the genome_id
						genome_id = host_spacer_pairs[a][i][0]
						genome_id = genome_id.split("|")[0]
						rest_id = host_spacer_pairs[a][k][0].split(genome_id)[1:]

						if (rest_id != []):
							host_spacer_pairs[a][k][0] = genome_id + rest_id[0]
						else:
							host_spacer_pairs[a][k][0] = genome_id
							host_spacer_pairs[a][k][-4] = host_spacer_pairs[a][i][-4]
					
					elif ((host_spacer_pairs[a][k][0], host_spacer_pairs[a][k][-4]) in array_names):
						array_names[(host_spacer_pairs[a][i][0],host_spacer_pairs[a][i][-4])] = (host_spacer_pairs[a][k][0], host_spacer_pairs[a][k][-4])
						genome_id = host_spacer_pairs[a][k][0]
						genome_id = genome_id.split("|")[0]
						rest_id = host_spacer_pairs[a][i][0].split(genome_id)[1:]
						if (rest_id != []):

							host_spacer_pairs[a][i][0] = genome_id + rest_id[0]
						else:
							host_spacer_pairs[a][i][0] = genome_id
						host_spacer_pairs[a][i][-4] = host_spacer_pairs[a][k][-4]
					else:
						array_names[(host_spacer_pairs[a][i][0],host_spacer_pairs[a][i][-4])] = (host_spacer_pairs[a][i][0], host_spacer_pairs[a][k][-4])
						genome_id = host_spacer_pairs[a][i][0]
						genome_id = genome_id.split("|")[0]
						rest_id = host_spacer_pairs[a][k][0].split(genome_id)[1:]

						if (rest_id != []):
							host_spacer_pairs[a][k][0] = genome_id + rest_id[0]
						else:
							host_spacer_pairs[a][k][0] = genome_id
						host_spacer_pairs[a][k][-4] = host_spacer_pairs[a][i][-4]

				k += 1
			i += 1
		a += 1
	for group in host_spacer_pairs:
		for row in group:
			# returns every element of the row
			ret_host_spacer_pairs.append(row)			
	return ret_host_spacer_pairs

def distance_based_deduplication(host_spacer_pairs,allowed_distance_between_pairs=50):
	ret_host_spacer_pairs = []	
	distance_matrix = []
	a = 0
	# iterate through each possible host target pairing
	while (a < len(host_spacer_pairs)):
		i = 0
		# allowed overlap len(spacer) - 2 * kmer
		# if the overlap is negative or 0 then this is allowed
		while (i < len(host_spacer_pairs[a])):
			
			allowed_overlap = 20
			if (float(host_spacer_pairs[a][0][-2]) < allowed_overlap):
				host_spacer_pairs[a].pop(0)
				continue
			k = i + 1
			while (k < len(host_spacer_pairs[a])):
				
				host_spacer_distance = abs(float(host_spacer_pairs[a][i][-2]) - float(host_spacer_pairs[a][k][-2])) # index should be the mapped distance.
				query_start = int(host_spacer_pairs[a][k][6])
				query_end = int(host_spacer_pairs[a][k][7])
				query_size = abs((query_start - query_end))
				query_length = int(host_spacer_pairs[a][k][5])
				# could change this to 20.
			#	allowed_overlap = query_length - 2 * query_size
			#	allowed_overlap = 20
				if (float(host_spacer_pairs[a][k][-2]) < allowed_overlap or host_spacer_distance < allowed_overlap):
					host_spacer_pairs[a].pop(k)

				else:
					k += 1
			i += 1
		a += 1	
	for group in host_spacer_pairs:
		for row in group:
			# Somehow returns every element of the row
			ret_host_spacer_pairs.append(row)			
	return ret_host_spacer_pairs

# deduplicate sequences by distance based indexing
# need to modify to take into account pairwise distance combinations
def distance_deduplication(host_spacer_pairs,allowed_distance_between_pairs=50):
	ret_host_spacer_pairs = []
	# iterate through each possible host target pairing.
	distance_matrix = []
	distance_row_index = -2
	a = 0
	while a < len(host_spacer_pairs):
		distance_dict = {}
		i = 0
		while i < len(host_spacer_pairs[a]):
			distance = abs(float(host_spacer_pairs[a][i][distance_row_index]))
			for dict_distance in distance_dict.keys():
				if (abs(distance -  dict_distance) < allowed_distance_between_pairs):
					break
			else:
				if (distance > allowed_distance_between_pairs):
					distance_dict[distance]	= host_spacer_pairs[a][i]
			i += 1
		for row in distance_dict.values():
			ret_host_spacer_pairs.append(row)	
		a += 1

	return ret_host_spacer_pairs

with open(sys.argv[1],"r") as csvfile:
	mapped_spacers = list(csv.reader(csvfile))


# will need to dictionalise by array or phage. The actual distance elimination could be done by a shared function
# first do dictionalisation by phage to eliminate the case of two arrays mapping to the same spacer.

phage_dict = {}
for spacer in mapped_spacers[1:]:
	if (spacer[1] not in phage_dict):
		phage_dict[spacer[1]] = [spacer]
	else:
		phage_dict[spacer[1]].append(spacer)	

# dictionalisation by phage to eliminate the case of two arrays mapping to the same spacer
mapped_spacers = distance_based_deduplication_array_names(list(phage_dict.values()))
print("mapped_spacers")

ret_out = open(sys.argv[2],"w")
spamwriter = csv.writer(ret_out)
for spacer in mapped_spacers:
	spamwriter.writerow(spacer)
ret_out.close()
array_dict = {}
for spacer in mapped_spacers[1:]:
	array_id = (spacer[0].split("|")[0],spacer[-4])
	if (array_id not in array_dict):
		array_dict[array_id] = [spacer]
	else:
		array_dict[array_id].append(spacer)
print("Doing distance deduplication now!")

# dictionalisation by array to eliminate multiple spacer mappings to redundant phage genomes.
mapped_spacers = distance_deduplication(list(array_dict.values()))

ret_out = open(sys.argv[3],"w")
spamwriter = csv.writer(ret_out)
for spacer in mapped_spacers:
	spamwriter.writerow(spacer)
ret_out.close()	











