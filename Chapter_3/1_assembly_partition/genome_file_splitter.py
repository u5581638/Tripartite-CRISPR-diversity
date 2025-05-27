from Bio import SeqIO
import sys
import os

# program to split concatenated genome assemblies into smaller blocks.
# INPUT: Concatenated ~10TB FASTA file of genome-assemblies
# OUTPUT: ~250-300mb files of FASTA sequences which are split from the input.
# i.e. genome_block_0 .. genome_block_1 .. genome_block_2 ..etc )
# SHELL: python3 genome_file_splitter.py {INPUT}

FILE = open(sys.argv[1], "r")
# FILE.seek(5008028911291)
sequence_iter = SeqIO.parse(FILE, "fasta")
i = 0
out_file_name = "genome_block_"

for sequence in sequence_iter:
	f = open(out_file_name + str(i) +".fasta", "a")
	if (os.path.getsize(out_file_name + str(i) + ".fasta") > 250000000): # get size > 250-500mb
	  	i += 1
	  	f.close()
	  	print(i)
	else:
		SeqIO.write(sequence, f, "fasta")