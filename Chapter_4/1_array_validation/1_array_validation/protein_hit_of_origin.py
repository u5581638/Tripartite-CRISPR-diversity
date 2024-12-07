from Bio import SeqIO
import csv
import sys

# retrieve protein sequences from hits to a BLASTp search result table (in csv format)
def anno_list_from_hits (orfs, hits):
	sequences = SeqIO.parse(orfs, "fasta")
	with open(hits, "r") as csvfile:
		hit_table = list(csv.reader(csvfile))
	out_file = open (orfs + "_annotated.fasta", "w")
	spam_writer = csv.writer(out_file)
	seq_dict = {}
	for sequ in sequences:
		seq_dict[sequ.id] = sequ  

	for hit in hit_table:
		if (hit[1] in seq_dict):
			SeqIO.write(seq_dict[hit[1]], out_file, "fasta")
	out_file.close()
	return orfs + "_annotated.fasta"