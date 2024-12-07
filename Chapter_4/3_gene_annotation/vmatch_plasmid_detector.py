# vmatch_plasmid_detector
# program to indicate possible vmatch matches

import sys
import csv
from Bio import SeqIO
# input maximal repeats from vmatch
# need the plasmid sequence to determine the length as well as the 
def line_denuller (split_sequences):
	ret_list = []
	for ele in split_sequences:
		if (ele != ''):
			ret_list.append(ele)
	return ret_list

# detect sequence repeats on contig indicative of circular DNA. Note: THis did not end up being used.
def plasmid_detect(vmatch_match_url_i, sequence_url, output_url):
	direct_repeat_threshold = 0
	plasmid_dr_threshold = 30
	distance_to_edge = 5

	sequence = SeqIO.read(sequence_url, "fasta")
	vmatch_output_url = open(vmatch_match_url_i, "r")

	v_match_output = vmatch_output_url.read()
	v_match_output = v_match_output.split("# find direct substring matches (repeats)\n")
	if (len(v_match_output)  > 1): # may simply not split. this will be a common source of bugs
		v_match_output = v_match_output[1]
	else:
		v_match_output = v_match_output[0]
	v_match_output = v_match_output.split("# overall space peak")
	v_match_output = v_match_output[0]
	v_match_output = v_match_output.split('\n')
	i = 0 
	while (i < len(v_match_output)):
		v_match_output[i] = line_denuller(v_match_output[i].split(" "))
		i += 1
#	print(v_match_output)
	
	i = 1 
	while (i < len(v_match_output)):
		if (len(v_match_output[i]) < 1):
			i += 1
			continue
	#	print(len(v_match_output))	
		if (int(v_match_output[i][2]) - distance_to_edge < 0 and  int(v_match_output[i][6]) + v_match_output[i][4] + distance_to_edge > len(sequence)): 
			print("target may be plasmid!") # instead may want to save the vmatch DNA segment to file
		i += 1

	i = 1	

	ret_out = open(output_url + "_drs" + ".csv", "w" )
	spam_writer = csv.writer(ret_out)
	spam_writer.writerow(["Genome_id", "length", "sequence_number", "relative_position", "type", "length2", "sequence_number2", "relative_position2","distance_value","E-value","score","perc-identity"])
	while (i < len(v_match_output)):
		if (len(v_match_output[i]) < 1):
			print(len(v_match_output[i]))
			i += 1
			continue

		if (len(v_match_output[i][0]) > direct_repeat_threshold or len(v_match_output[i][4]) > direct_repeat_threshold):
			spam_writer.writerow([sequence.id] + v_match_output[i])
		i += 1

	ret_out.close()
	vmatch_output_url.close()
	return 0