# mapped_spacer_entry_matching
# script to match kmers with their corresponding mapped spacer entries

import csv
import sys

with open(sys.argv[1],"r") as csvfile:
	spacer_mapping_table = list(csv.reader(csvfile))

with open(sys.argv[2],"r") as csvfile:
	detected_kmers = list(csv.reader(csvfile))
# create a dictionary of mapped_spacers
spacer_mapping_dict = {}
for spacer in spacer_mapping_table[1:]:
	if ((spacer[0].split("|")[0], spacer[-2]) not in spacer_mapping_dict):
		spacer_mapping_dict[(spacer[0].split("|")[0], spacer[-2])] = [spacer]
	else:
		spacer_mapping_dict[(spacer[0].split("|")[0], spacer[-2])].append(spacer)	


kmer_dict = {}
for kmer in detected_kmers[1:]:
	if ((kmer[0].split("|")[0],kmer[-2]) not in kmer_dict):
		kmer_dict[(kmer[0].split("|")[0],kmer[-2])] = [kmer]

	else:
		kmer_dict[(kmer[0].split("|")[0],kmer[-2])].append(kmer)
for kmer in kmer_dict:
	if (kmer in spacer_mapping_dict):
		kmer_dict[kmer] = kmer_dict[kmer] + spacer_mapping_dict[kmer]
		print("joined")
		# Transfer and transform the information on each whole spacer into coordinates on the mapped phage contig

ret_out = open(sys.argv[3],"w")
spamwriter = csv.writer(ret_out)

for kmer in kmer_dict:
	for row in kmer_dict[kmer]:

		spamwriter.writerow(row)
ret_out.close()