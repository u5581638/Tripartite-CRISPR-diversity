import csv
import sys
from Bio import SeqIO

# script to retrieve sequence metadata for JGI-derived sequences using a merged metadata table downloaded from the GOLD database website
# SHELL: python3 crispr_window_metadata_retrieval.py true_gold_compile3_merged.csv spliced_debug_corrected_cleaned_rerun_20kb_local_windows_combined.fasta output file name>
# INPUT: GOLD table (merged)
# i.e. true_gold_compile3_merged.csv
with open (sys.argv[1],"r") as csvfile:
	microbiome_table = list(csv.reader(csvfile))

micro_biome_dict = {}

for row in microbiome_table[1:]:
	micro_biome_dict[row[12]] = row 

# DNA within 20kb of a CRISPR-array
# i.e. spliced_debug_corrected_cleaned_rerun_20kb_local_windows_combined.fasta
sequences = SeqIO.parse(sys.argv[2],"fasta")

# OUTPUT: GOLD annotation data for each contig with a Ga number in each 40kb sequence window
ret_out = open(sys.argv[3],"w")
spamwriter = csv.writer(ret_out)

for sequ in sequences:
	sequ_id = sequ.id.split("::")[0]
	sequ_id = sequ.id.split("_")[0]
	if sequ_id in micro_biome_dict:
		spamwriter.writerow(micro_biome_dict[sequ_id])

ret_out.close()
