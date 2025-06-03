from Bio import SeqIO
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Get Gaetan to review this code!!

# script to deduplicate raw contigs by 
# Use the cleaned contigs file present on the NCI as input
# This script will need to run on Gaetan's computer.

# may want to expand upon this function to return values which include the start and end coordinates of the genomes!!!!
def tuple_sorter (contig_tuple):
	return int(contig_tuple[1])

def window_set_generator (contigs):
	window_set_dict = {}
	for contig in contigs:
		contig_id = contig.id.split("::")
		contig_index = contig_id[1]
		contig_id = contig_id[0]
		contig_index = contig_index.split(":")
		contig_start_index = contig_index[0]
		contig_end_index = contig_index[1]
		# here write code to extract the start and end coordinates

		if (contig_id not in window_set_dict):
			window_set_dict[contig_id] = [[contig, contig_start_index, contig_end_index]]
		else:
			window_set_dict[contig_id].append([contig, contig_start_index, contig_end_index])
	return window_set_dict	

def amalgamate (previous_contig, contig):
	#need to splice the exta sequence into the contig.
	starting_splice_index_next_contig = int(previous_contig[2]) - int(contig[1])
	end_splice_index_next_contig = int(contig[2])
	splice_contig = str(contig[0].seq)
	splice_contig[starting_splice_index_next_contig:end_splice_index_next_contig]
	new_contig_seq = Seq(str(previous_contig[0].seq) + splice_contig)
	new_contig = contig[0]
	new_contig_id = new_contig.id.split("::")
	new_contig_id = new_contig_id[0]
	new_contig_id = new_contig_id + "::" + str(previous_contig[1]) + ":" + str(end_splice_index_next_contig)
	new_contig.id = new_contig_id
	new_contig.seq = new_contig_seq
#	print ("The original contig sequence is:")
#	print (previous_contig[0].seq)
#	print ("The 2nd original contig sequence is:")
#	print (contig[0].seq)
#	print("The joined contig sequence is:")
#	print (new_contig)
	new_contig = [new_contig, previous_contig[1], contig[2]]
	return new_contig	

def longest_overlapping_contig (contigs, file_handle):

	ret_out = open(file_handle, "a")
	# sort contigs by start index
	contigs.sort(key=tuple_sorter)
	previous_contig=contigs[0]
	print(len(contigs))
	i = 1
	while(i < len(contigs)):
		if ( int(previous_contig[2]) >= int(contigs[i][1]) ):
			previous_contig = amalgamate(previous_contig, contigs[i])
		else: # contig not 
			print("Hello")
			SeqIO.write(previous_contig[0], ret_out, "fasta")
			previous_contig = contigs[i] 
		i += 1	
	print("Done!!")	
	print(previous_contig)	
	SeqIO.write(previous_contig[0], ret_out, "fasta")
	ret_out.close()
	return 0	

	# check whether the first contig overlaps with the next contig and make the next contig the target contig
	# if yes
	# amlagamate
	# if no
	# make the next contig the targt contig
	# write all sequences to memory/file. This will consume immense amounts of memory/IO. Sequences can be written whenever the no statement is called.

contigs = SeqIO.parse(sys.argv[1], "fasta")

window_set_dict = {} # sort the contigs into a list of individual windows. The key should be the genome id. The values should be the windows from each genome.
window_set_dict = window_set_generator(contigs)

windows = list(window_set_dict.values())
#print(windows)
print("2nd_window:")
print(len(windows))
# need to go through the windows for each genome and extract the minimum and maximum values:
# first need to check if any of the windows overlap.
# most efficent way to do this?

for window in windows: 
	longest_overlapping_contig(window, "spliced_" + sys.argv[1])

