# protein_retriever.py

from Bio import SeqIO
import sys
import csv

# INPUT: 1. results from BLASTp (as csv table)
# 		 2. FASTA file corresponding to the BLASTp database used
# OUTPUT: retrieved protein sequences (in FASTA format)

proteins = SeqIO.parse(sys.argv[2],"fasta")
prot_dict = {}

for prot in proteins:
	prot_dict[prot.id] = prot

with open (sys.argv[1], "r") as csvfile:
	hit_table = list(csv.reader(csvfile))

ret_list = []

for hit in hit_table:
	if (hit[1] in prot_dict):
		ret_list.append(prot_dict[hit[1]])
SeqIO.write(ret_list,sys.argv[3],"fasta")	