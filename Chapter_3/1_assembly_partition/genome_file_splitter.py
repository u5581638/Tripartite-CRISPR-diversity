from Bio import SeqIO
import sys
import os

FILE = open(sys.argv[1], "r")
FILE.seek(5008028911291)
sequence_iter = SeqIO.parse(FILE, "fasta")
i = 20032
out_file_name = "genome_block_"

for sequence in sequence_iter:
	f = open(out_file_name + str(i) +".fasta", "a")
	if (os.path.getsize(out_file_name + str(i) + ".fasta") > 250000000): # get size > 250-500mb
	  	i += 1
	  	f.close()
	  	print(i)
	else:
		SeqIO.write(sequence, f, "fasta")