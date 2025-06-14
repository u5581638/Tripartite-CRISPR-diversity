import csv
import sys
from Bio import SeqIO
import re

# script to retrieve sequence metadata for NCBI-derived sequences using a merged metadata table downloaded from the GOLD database website
# SHELL: python3 crispr_window_metadata_retrieval_ncbi.py master_genbank_table.csv 
# INPUT: master genbank table from NCBI FTP repository with additional entries added by entrez queries.
# i.e. master_genbank_table.csv spliced_debug_corrected_cleaned_rerun_20kb_local_windows_combined.fasta <output file name>
with open (sys.argv[1],"r") as csvfile:
	microbiome_table = list(csv.reader(csvfile))

micro_biome_dict = {}

for row in microbiome_table[1:]:
#	print(row)
	id_row = row[3] #need to split on the characters.
	id_row = re.match('[A-Z]*',id_row)
	id_row = id_row.group(0)
	micro_biome_dict[id_row] = row 
# INPUT: sequences with NCBI identifers from sequences within 20kb of a CRISPR-array
# spliced_debug_corrected_cleaned_rerun_20kb_local_windows_combined.fasta
sequences = SeqIO.parse(sys.argv[2],"fasta")
# OUTPUT: NCBI annotation data for each contig with a matching NCBI ID in each 40kb sequence window
ret_out = open(sys.argv[3],"w")
spamwriter = csv.writer(ret_out)



for sequ in sequences:
	identifier = sequ.id.split("::")[0]
	if (re.search('.*\.[0-9]',identifier) and not re.search('Contig.*$',identifier) and not re.search('NODE.*$',identifier) and not re.search('_',identifier)):
		ncbi_id = re.match('[A-Z]*',identifier)
		ncbi_id = ncbi_id.group(0)
		if ncbi_id in micro_biome_dict:
			spamwriter.writerow(micro_biome_dict[ncbi_id])

ret_out.close()
