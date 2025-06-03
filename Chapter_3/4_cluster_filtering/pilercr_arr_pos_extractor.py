# script to extract the pilercr id and start/end sites in a usable form for filtration

from Bio import SeqIO 
import sys
import re
import csv

def line_denuller (split_sequences):
	ret_list = []
	for ele in split_sequences:
		if (ele != ''):
			ret_list.append(ele)
	return ret_list

# INPUT: file containing CRISPR-array predictions in quasi-FASTA format
# i.e. spliced_cleaned_rerun_20kb_local_windows_combined.fasta_pilercr.fa 
# OUTPUT: table containing with each row giving the contig_id and array start and end position.
# i.e. spliced_cleaned_rerun_20kb_local_windows_combined.fasta_pilercr.fa_real_arr_positions.csv
# SHELL: python3 pilercr_arr_pos_extractor.py spliced_cleaned_rerun_20kb_local_windows_combined.fasta_pilercr.fa

characters = re.compile('[0-9]') 
equals = re.compile('=')
repeat_pattern = re.compile('') 
sequence_file = open(sys.argv[1], "r") 
sequence_mega_str = sequence_file.read()
sequences = sequence_mega_str.split(">") 
ret_file = open(sys.argv[1] + "_real_arr_positions.csv", "a")
switch = False
#print(sequences)
sequences = sequences[1:]
for sequence in sequences:
#	print(sequence)
	sequence_lines = sequence.split("\n")
	line = sequence_lines[0]
	sequence_lines = sequence_lines[1:]
	header = line # first line
	ret_lines = []
	#print(sequence_lines)
	for line in sequence_lines:
		line_items = line.split(" ")
		line_items = line_denuller(line_items)
		ret_line = []
		if (line_items != []):
			print(line_items[0][0])
		if (line_items == []):
			continue
			
		elif (switch == False and re.search(equals , line_items[0][0])): 

			switch = True
			continue
		elif (switch == True and re.search(characters, line_items[0][0])):

			ret_line.append(line_items[0])
			ret_line.append(line_items[1])
			ret_lines.append(ret_line)
			continue
		elif (switch == True and re.search(equals, line_items[0][0])):
			switch = False			
			break
		else:
			continue	
	start_pos = ret_lines[0][0]
	end_pos = str(int(ret_lines[-1][1]) + int(ret_lines[-1][0]))
	spam_writer = csv.writer(ret_file)
	spam_writer.writerow(["".join(header)] + [start_pos] + [end_pos])
ret_file.close()



