import sys
from Bio import SeqIO
import re
import copy


 
# file containing genome_aseemblies in FASTA format
fasta_sequences = list(SeqIO.parse(sys.argv[1], "fasta"))
# file containing PILER-CR crispr array predictions in FASTA format.
crispr_arrays = list(SeqIO.parse(sys.argv[2], "fasta"))
window_size = sys.argv[3]

# convert to dict 
# 1. go through all the sequences
# 2. go through all the crispr arrays
# 3. run extraction function
# write extracted windows to file
def window_extractor (genome, index, window_size):
	genome_copy = genome
	if(index - window_size < 0 and index + window_size > len(genome_copy)):
		genome_copy.id = genome_copy.id + "::0:" +str(len(genome_copy))
		return genome_copy
	elif (index - window_size < 0):
		genome_copy.id = genome_copy.id + "::0:" + str((index+window_size))
		return genome_copy[:(index+window_size)]
	elif (index + window_size > len(genome_copy)):
		genome_copy.id = genome_copy.id + "::" + str((index-window_size)) + ":" + str(len(genome_copy))
		return genome_copy[(index-window_size):]
	else:
		genome_copy.id = genome_copy.id + "::" + str((index-window_size)) + ":" + str((index+window_size))
		return genome_copy[(index-window_size):(index+window_size)] 

i = 0
position_pattern = re.compile("Pos=.*]")
#print (len(fasta_sequences))
#print(len(crispr_arrays))
extracted_windows = []
#j = 0
while (i < len(fasta_sequences)):
	k = 0
	while (k < len(crispr_arrays)):
	#	print(k)
		crispr_arrays_id = crispr_arrays[k].description.split("[")
		crispr_arrays_id = crispr_arrays_id[0]
		crispr_arrays_id = crispr_arrays_id.split(" ")
		crispr_arrays_id = crispr_arrays_id[0]

		if (fasta_sequences[i].id == crispr_arrays_id):
			matched_str = re.search(position_pattern, crispr_arrays[k].description)
			matched_str = str(matched_str.group(0))
			crispr_index = int(matched_str[4:-1])
			sequences_id_copy = copy.deepcopy(fasta_sequences[i])
		#	print(sequences_id_copy.id)
		#	print(fasta_sequences[i].id)
			sequence_window = window_extractor(sequences_id_copy, crispr_index, int(window_size))
		#	print(fasta_sequences[i].id)
		#	print(j)
			extracted_windows.append(sequence_window)
		#	print(extracted_windows)
#			j += 1	
		k += 1
	i += 1 
SeqIO.write(extracted_windows, sys.argv[1] + "_windows.fasta", "fasta")