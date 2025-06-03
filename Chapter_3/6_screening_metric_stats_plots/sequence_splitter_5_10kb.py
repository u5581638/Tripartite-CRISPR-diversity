import sys
from Bio import SeqIO

# program to split FASTA file of representative protein sequences into individual sequences.
# INPUT: FASTA file of representative protein sequences (5-10kb)
# i.e. all_run_5_10kb_sequences.fasta
# OUTPUT: individual sequences in FASTA file <- labelled by id
# i.e. "split_sequences_5_10kb/*"
# SHELL: python3 sequence_splitter_5_10kb.py all_run_5_10kb_sequences.fasta
sequences = SeqIO.parse(sys.argv[1],"fasta")

for sequ in sequences:
	seq_file = sequ.id.split("|")
	seq_file = "sequence_tmp_" + seq_file[0] + seq_file[1]
	SeqIO.write(sequ, "split_sequences_5_10kb/" + seq_file,"fasta")