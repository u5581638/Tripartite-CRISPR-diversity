# host_header_adder

# script to add the host genome to the header names in each FASTA sequence

import csv
import sys
import subprocess
from Bio import SeqIO
import os
# need to first retrieve genomes

with open(sys.argv[1], "r") as csvfile:
	host_phage_hitmap = list(csv.reader(csvfile))



# break into individual host_genome groups
genome_host_dict = {}
genome_ids = set()
phage_ids = set()
b = "cas12g.fasta_all_hits.csv_genomes.fasta"
b_phage = "cas12g.fasta_all_hits.csv_genomes.fasta_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_hitmap.csv_standardised.csv_non_redundant.csv_filtered.csv_deduplicated.csv_filtered_hits_extracted_faidx_bp_window.fasta"
out_file_name = "cas12g_retrieved_genomes.fasta"
out_file_name2 = "cas12g_retrieved_genomes_phage.fasta"
out_file = open(out_file_name2,"a")
for row in host_phage_hitmap[1:]:
	if (row[0] not in genome_ids):
		genome_ids.add(row[0])
		subprocess.run(["samtools faidx " + b + " " + "\'" + row[0] + "\'" + " >> " + out_file_name], shell=True)
		
	if (row[4] not in phage_ids):
		phage_ids.add(row[4])
		subprocess.run(["samtools faidx " + b_phage + " " + "\'" + row[4] + "\'" + " > " + "phage_file.fa"], shell=True)
		sequence = SeqIO.read("phage_file.fa","fasta")
		sequence.id = sequence.id + "|" + "host_genome=" + row[0]
		sequence.description = ""
		SeqIO.write(sequence,out_file,"fasta")

out_file.close()		

# subprocess.run(["prodigal -i " + out_file_name + " -o " + out_file_name + ".txt" + " -a " + out_file_name + "_aa_raw.fasta"],shell=True)
# subprocess.run(["prodigal -i " + out_file_name2 + " -o " + out_file_name2 + ".txt" + " -a " + out_file_name2 + "_aa_raw.fasta"],shell=True)