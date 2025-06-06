# phage_contig_start_adder.py

# add phage contig start to kmer_entries.

# need to change the distance to the minimum number during spacer_expansion

import sys
import csv
from Bio import SeqIO
import re

# INPUT: 1. mapped target site contigs (phages) in FASTA file.
#		 i.e. cas12a.fasta_all_hits.csv_genomes.fasta_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_hitmap.csv_standardised.csv_non_redundant.csv_filtered.csv_deduplicated.csv_filtered_hits_extracted_faidx_bp_window.fasta
#		 2. file containing kmer matches to contigs (self-array targets removed)
#		 i.e. cas12a_detection_parallel_gadi_kmer14_29-3-2024.csv_non-self.csv
# OUTPUT: file containing kmer matches to contigs with contig start coordinate appended to the headers.
#		 i.e. cas12a_detection_parallel_gadi_kmer14_29-3-2024.csv_non-self.csv_pcontig.csv
# SHELL: python3 phage_contig_start_adder.py cas12a.fasta_all_hits.csv_genomes.fasta_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_hitmap.csv_standardised.csv_non_redundant.csv_filtered.csv_deduplicated.csv_filtered_hits_extracted_faidx_bp_window.fasta cas12a_detection_parallel_gadi_kmer14_29-3-2024.csv_non-self.csv cas12a_detection_parallel_gadi_kmer14_29-3-2024.csv_non-self.csv_pcontig.csv
phage_contigs = SeqIO.parse(sys.argv[1],"fasta")


phage_dict = {} # This only needs to carry the phage_contig_start length as the value. Use the existing legnth value in the table to compute the end of contig length
# need to overwrite one of the pre-existing table rows to store this information without altering the table dimensions
for contig in phage_contigs:
#	print(contig.id)
	if (len(contig.id.split(":")) < 2):
		continue
	contig_id = contig.id.split(":") [1]
	contig_start_site = contig_id.split("-") [0]
	contig_id = contig.id.split(":")[0]
	if (contig_id) not in phage_dict:
		phage_dict[contig_id] = [contig]
	else:
		phage_dict[contig_id].append(contig)	

print(phage_dict.keys())
print("DONE")
kmer_handle = open(sys.argv[2],"r")
kmer_table = csv.reader(kmer_handle)
header = next(kmer_table)
ret_out = open(sys.argv[3],"a")
spamwriter = csv.writer(ret_out)
spamwriter.writerow(header)

for kmer in kmer_table:
	switch = 0 
	target_start = kmer[8]
	target_end = kmer[9]
	phage_id = kmer[1]
	kmer_frag = kmer[4]
	# iterate over the possible phage fragments
	if (phage_id not in phage_dict):
		print("phage_id")
		print(phage_id)
		continue
	print("enter loop")
	for phage in phage_dict[phage_id]:
		phage_start_coord = int(phage.id.split(":") [1].split("-") [0])
		max_phage_len = len(phage) + phage_start_coord 
		# 
		if (int(phage_start_coord) < int(target_start) < int(max_phage_len) and int(phage_start_coord) < int(target_end) < int(max_phage_len)):
			my_contig = str(phage.seq).upper()
			real_target_start = int(target_start) - phage_start_coord
			real_target_end = int(target_end) - phage_start_coord
			if (real_target_start < real_target_end):
				# Should only scan each kmer mapped region. This should significantly increase the speed of the search
				m = re.search(phage[4],my_contig[real_target_start:real_target_end])
			else:
				m = re.search(phage[4],my_contig[real_target_end:real_target_start])
			if (m):
				spamwriter.writerow(kmer[:10] + [phage_start_coord] + [str(real_target_start)] + [str(real_target_end)] + kmer[13:])
				switch = 1
				break
	if (switch != 1):
		print("Mistake has occurred.")

ret_out.close()
kmer_handle.close()