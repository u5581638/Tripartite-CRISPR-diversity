# sequence_name_adder

# script to add the sequence name to the sequence files.

# INPUT: 5_10kb protein sequence annotations from pfam/padloc. 0_5kb uses the difference type of sequence_id with a different format.
# i.e. sequence_tmp_genome_block_9998.fasta_windows.fastagene_1710499_17_pfam_descriptions.txt_pfam_descriptions.txthits.csv
# OUTPUT: File with annotation added
# i.e. sequence_tmp_genome_block_9998.fasta_windows.fastagene_1710499_17_pfam_descriptions.txt_pfam_descriptions.txthits.csvheader_added.csv
# SHELL: python3 sequence_name_adder.py sequence_tmp_genome_block_9998.fasta_windows.fastagene_1710499_17_pfam_descriptions.txt_pfam_descriptions.txthits.csv

import csv
import sys

with open(sys.argv[1], "r") as csvfile:
	hit_table = list(csv.reader(csvfile))

ret_out = open(sys.argv[1] + "header_added.csv","w")
spamwriter = csv.writer(ret_out)
header = sys.argv[1].split("/")[1]
header = header.split("_")
header = header[0] + "_" + header[1] + "_" + header[2] + "_" + header[3] + "_" + header[4] + "_" + header[5] + "_" + header[6] + "_" + header[7]
for row in hit_table:
	out_row = [header] + row
	spamwriter.writerow(out_row)
ret_out.close()