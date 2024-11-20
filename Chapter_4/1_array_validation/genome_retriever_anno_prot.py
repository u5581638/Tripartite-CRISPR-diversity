from Bio import SeqIO
import sys
import csv


def genome_retriever (seq_url, hit_url):
	sequences = SeqIO.parse(seq_url, "fasta")
	with open(hit_url, "r") as csvfile:
		hit_table = list(csv.reader(csvfile))
	seq_dict = {}
	for sequ in sequences:
		seq_dict[sequ.id] = sequ
	hits = set(list(map(lambda x: x[1], hit_table)))
	ret_list = []
	for hit in hits:
		if hit in seq_dict:
			ret_list.append(seq_dict[hit])
		else:
			print("Error!! ID not recognised!!")
	return ret_list


