import sys
from Bio import SeqIO

# program to split FASTA file of representative protein sequences into individual sequences.
# INPUT: FASTA file of representative protein sequences (0-5kb)
# i.e. all_run_0_5kb_sequences.fasta
# OUTPUT: individual sequences in FASTA file <- labelled by id
# i.e. "split_sequences_0_5kb/*"
# SHELL: python3 sequence_splitter_0_5kb.py all_run_0_5kb_sequences.fasta

sequences = SeqIO.parse(sys.argv[1],"fasta")

for sequ in sequences:
	seq_file = sequ.id.split("|")
	seq_file = seq_file[0]
	SeqIO.write(sequ, "split_sequences_0_5kb/" + seq_file,"fasta")