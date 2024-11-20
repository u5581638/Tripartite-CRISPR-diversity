# program to reformat phage genome sequences

from Bio import SeqIO
import sys
import re
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord

def line_appender (header_str):
	ret_str = header_str[0]
	i = 1
	while (i < len(header_str) - 1):
		ret_str = ret_str + "-" + header_str[i]
		i += 1
	return ret_str

def reformat(spacer_url, output_dir):
	sequences = list(SeqIO.parse(spacer_url, "fasta"))
#	print(sequences)
	ret_list = []
	i = 0
	while (i < len(sequences)):
		header = sequences[i].id.split("-")
	#	print(header)
		end_sequence = header[1]

		if (int(end_sequence) > len(sequences[i])):

			sequences[i].id = line_appender(header) + '-' + end_sequence
			print(sequences[i].id)
		ret_list.append(sequences[i])
		i += 1
	SeqIO.write(ret_list, output_dir, "fasta")
	return 0

