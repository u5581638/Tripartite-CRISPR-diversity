# protein retriever

import csv
import sys
from Bio import SeqIO


sequences = SeqIO.parse(sys.argv[1],"fasta")
seq_dict = {}
for sequ in sequences:
	if (sequ.id not in seq_dict):
		seq_dict[sequ.id] = sequ

with open(sys.argv[2],"r") as csvfile:
	hit_table = list(csv.reader(csvfile))

row_set = set()
for row in hit_table:
	if row[1] not in row_set:
		row_set.add(row[1])

ret_list = []
for my_id in row_set:
	ret_list.append(seq_dict[my_id])
SeqIO.write(ret_list,sys.argv[3],"fasta")
