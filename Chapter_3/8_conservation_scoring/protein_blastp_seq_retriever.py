# retrieve sequences from 10TB block
import sys
import csv
from Bio import SeqIO
import fcntl

# INPUT: 1. Protein hits from query translated genome blocks in labelled_proteins/
#  		 2. Protein_translated genome blocks

# Output: 1. Protein sequences retrieved as a result of BLASTp search which are to a specific HMM based protein annotation
# SHELL: find /g/data/va71/labelled_proteins/ -name "*_aa.fa" -type "f" -exec basename {} \; | xargs -n 1 -I {} -P 104 python3 protein_blastp_seq_retriever.py uma2/{}_hits.csv /g/data/va71/labelled_proteins/{} uma2_proteins/{}_proteins.fasta

seq_hits_handle = open(sys.argv[1],"r")
hit_table = csv.reader(seq_hits_handle)
genome_block = SeqIO.parse(sys.argv[2],"fasta")
genome_dict = {} 
for contig in genome_block:
	if (contig.id not in genome_dict):
		genome_dict[contig.id] = contig

ret_out = open(sys.argv[3],"a")

for row in hit_table:
	if (row[1] in genome_dict):
		SeqIO.write(genome_dict[row[1]],ret_out,"fasta")

ret_out.close()