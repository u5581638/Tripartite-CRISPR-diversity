# genome_block_relabeller.py
# script to label the sequences of an input block with the block nanmes

from Bio import SeqIO
import sys

sequences = SeqIO.parse(sys.argv[1], "fasta")
ret_seq = []
for sequence in sequences:
        this_sequence = sequence
        this_sequence.description = this_sequence.description + "|" + sys.argv[1]
        ret_seq.append(this_sequence)
SeqIO.write(ret_seq, sys.argv[1] + "_labelled.fasta", "fasta")