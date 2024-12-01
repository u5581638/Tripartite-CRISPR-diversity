import csv
import sys
from Bio import SeqIO

# script to retrieval sequence metadata for JGI-derived sequences using a merged metadata table downloaded from the GOLD database website

# INPUT: GOLD table (merged)
with open (sys.argv[1],"r") as csvfile:
	microbiome_table = list(csv.reader(csvfile))

micro_biome_dict = {}

for row in microbiome_table[1:]:
	micro_biome_dict[row[12]] = row 

# DNA within 20kb of a CRISPR-array
sequences = SeqIO.parse(sys.argv[2],"fasta")

ret_out = open(sys.argv[3],"w")
spamwriter = csv.writer(ret_out)

for sequ in sequences:
	sequ_id = sequ.id.split("::")[0]
	sequ_id = sequ.id.split("_")[0]
	if sequ_id in micro_biome_dict:
		spamwriter.writerow(micro_biome_dict[sequ_id])

ret_out.close()
