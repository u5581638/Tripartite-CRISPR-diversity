# algorithm to parse prodigal protein orf output files into genemark file format
# Seperate headers using |
# Order of headers:
# 1. Protein_orf_id
# 2. ORF prediction program (prodigal)
# 3. protein size
# 4. sense strand
# 5. start of ORF (relative to contig)
# 6. end of ORF (relative to contig)
# 7. genome_id
# 8. Organism name
# 9. sequence

from Bio import SeqIO
import sys
import csv

def sense (value):
	print (value)
	if (value == "-1"):
		return "-"
	elif (value == "1"):
		return "+"
	else:
		print ("error, wrong string parsed")
		return "error"		

# INPUT: list of protein sequence clusters sorted in descending order.
# i.e. sorted_descending_deduplicated_flattened_all_sequences.fasta
# OUTPUT: list of protein sequence clusters sorted in descending order, re-formatted
# i.e. sorted_descending_deduplicated_flattened_all_sequences.fasta_p_formatted.fasta
# SHELL: python3 prodigal_parser_cluster.py sorted_descending_deduplicated_flattened_all_sequences.fasta

sequences = list(SeqIO.parse(sys.argv[1] , "fasta"))

i = 0
while (i < len(sequences)):
	
	seq_description = sequences[i].description
	cluster_number = seq_description.split("|cluster")
	cluster_number = cluster_number[-1]
	seq_description_v1 = seq_description.split(" # ")
	seq_description_v2 = seq_description.split(";")
	seq_prot_id = seq_description_v2[0]
	seq_prot_id = seq_prot_id.split("=")
	print (seq_prot_id)
	seq_prot_id = seq_prot_id[1]
	sequences[i].id = "gene_" + seq_prot_id + "|" + "prodigal" + "|" + str(len(sequences[i].seq)) + "aa" + "|" + sense(seq_description_v1[3]) + "|" + seq_description_v1[1] + "|" + seq_description_v1[2] + "|" + seq_description_v1[0] + "|" + "cluster_" + cluster_number
	sequences[i].description = "gene_" + seq_prot_id + "|" + "prodigal" + "|" + str(len(sequences[i].seq)) + "aa" + "|" + sense(seq_description_v1[3]) + "|" + seq_description_v1[1] + "|" + seq_description_v1[2] + "|" + seq_description_v1[0] + "|" + "unknown_bacterial genome" + "|" + "cluster_" + cluster_number
	i += 1
SeqIO.write(sequences, sys.argv[1] + "_p_formatted.fasta", "fasta")	
