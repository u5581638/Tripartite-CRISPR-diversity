# program to deduplicate protein translations of contigs which have been inferred to be the source of all headaches

# do this the old fashioned way.
# if a given contig matches the sequence id, contig + aa positions (effective genome position) and sense direction, then it MUST be identical.
# first need to look at the structure of a prodigal output file.

from Bio import SeqIO
import sys

# INPUT: Predicted protein translations from prodigal in FASTA format, prior to sequence clustering using mmseqs2.
# i.e. all_sequences.fasta
# OUTPUT: Predicted proteins without any duplicate entries from duplicate file headers or duplicate sequence windows
# i.e. 1. deduplicated_all_sequences.fasta
#	   2. problematic_all_sequences.fasta (anomalous contig_id structure <- error catch case)
# SHELL: python3 contig_deduplicator.py all_sequences.fasta

sequences = SeqIO.parse(sys.argv[1], "fasta")
seq_dict = {} # want a list of sequences with identical identifers

problematic_sequences = []
for sequ in sequences:
	sequence_decomp = sequ.description.split(" # ") # keep .id actually works. May need to use .description
#	se
	sequence_id = sequence_decomp[0]
	sequence_id = sequence_id.split("::")
	genome_label = sequence_id[0]
#	print(sequ.id)
#	print (sequence_id)
	if (len(sequence_id) < 2):
		problematic_sequences.append(sequ)
	else:	
		sequence_id = sequence_id[1]
		sequence_id = sequence_id.split("_")
		sequence_id = sequence_id[0]
		sequence_id = sequence_id.split(":")
		start_contig_pos = sequence_id[0]
		end_contig_pos = sequence_id[1]

	# index for negative strands sequences on the starting coordinate because some contigs in all sequence data may not have an end coordinate.
	#print(sequence_decomp)
	if (sequence_decomp[3] == "1"):
		pos = str(int(start_contig_pos) + int(sequence_decomp[1]))
	else: # secquence_decomp[3] == "-1" 
		pos = str(int(end_contig_pos) + int(sequence_decomp[2]))
	sequence_id = genome_label + pos 
		

	if (sequence_id) not in seq_dict:
		seq_dict[sequence_id] = sequ
#		print("Gotcha!")
SeqIO.write(seq_dict.values(), "deduplicated_" + sys.argv[1], "fasta")		
SeqIO.write(problematic_sequences,  "problematic_" + sys.argv[1], "fasta")	
			

