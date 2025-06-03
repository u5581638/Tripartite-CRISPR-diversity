# program to parse the output of hmmscan for each result line

from Bio import SeqIO
import sys
import csv
import copy
import re

def line_denuller (split_sequences):
	ret_list = []
	for ele in split_sequences:
		if (ele != ''):
			ret_list.append(ele)
	return ret_list

def hmmscan_parse(file_name):

	in_file = open(file_name, "r")
	ret_list = []
	seq_hits = in_file.read() # type-casting to str should be unessessary!
	if re.search('>> ', seq_hits):
		seq_hits = seq_hits.split(">> ")
		seq_hits = seq_hits[1:]
	else:
		return []	
	outfile = open(file_name +"_hmmscan.csv", "w")
	spam_writer = csv.writer(outfile)
	ret_list = []
	equal3_switch = 0

	for sequ in seq_hits:
		all_lines = sequ.split("\n")
		info_line = all_lines[3]
		description = all_lines[0]
		my_id = description.split(" ")
		my_id = my_id[0]
		
		info_line = line_denuller(info_line.split(" "))
	#	print(my_id)
	#	print(info_line)
		score = info_line[2]
		c_evalue = info_line[4]
		i_evalue = info_line[5]
		hmm_start = info_line[6]
		hmm_end = info_line[7]
		ali_start = info_line[9]
		ali_end = info_line[10]
		ret_list.append([my_id, description, score, c_evalue, i_evalue, hmm_start, hmm_end, ali_start, ali_end])
		spam_writer.writerow([my_id, description, score, c_evalue, i_evalue, hmm_start, hmm_end, ali_start, ali_end])
	outfile.close()	
	return ret_list

# INPUT: Output of DEFLOC queried file from HHblits.
#  i.e. sequence_999_crisprcas_hits.txt
# OUTPUT: Conversion of output to csv format.
# i.e. sequence_997_crisprcas_hits.txt_hmmscan.csv
# SHELL: See running_scripts/pfam_annotations_running_script.sh

hmmscan_parse(sys.argv[1])

