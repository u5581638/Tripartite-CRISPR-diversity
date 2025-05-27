# genome_block_relabeller.py
# script to label the sequences of an input block with the block nanmes
# INPUT: One of 38810 assembled genome sequence data files in FASTA format. 
# i.e. labelled_genomes/genome_block_33.fasta
# OUTPUT: assembled genome sequence data file with headers including file name
# i.e. labelled_genomes/genome_block_33.fasta_labelled.fasta
# SHELL: find labelled_genomes/ -name "*" | python3 genome_block_relabeller.py {INPUT} 

from Bio import SeqIO
import sys

sequences = SeqIO.parse(sys.argv[1], "fasta")
ret_seq = []
for sequence in sequences:
        this_sequence = sequence
        this_sequence.description = this_sequence.description + "|" + sys.argv[1]
        ret_seq.append(this_sequence)
SeqIO.write(ret_seq, sys.argv[1] + "_labelled.fasta", "fasta")