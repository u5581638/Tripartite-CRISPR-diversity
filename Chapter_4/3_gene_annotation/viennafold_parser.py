# program to parse viennafold diagram output!
# Output table should show the free energy for each hairpin! These can be used to assess whether or not secondary structure exists!!

import sys
import csv

def line_denuller(input_str):
	ret_str = input_str[0]
	i = 1
	while (i < len(input_str) - 1):
		ret_str = ret_str + "_" + input_str[i]
		i += 1
	return ret_str


# parse RNA secondary structure predictions to table.
def parse(cds_url, ret_url):
	my_sequence_handle = open(cds_url, "r")
	my_sequence = my_sequence_handle.read()
	my_sequence_handle.close()
	my_sequence = my_sequence.split("\n") # split into each line
	i = 1
	header = my_sequence[0].split(" ") [0] # take just the protein id - check this!!!
	ret_list = []
	print(my_sequence)
	mea_line = my_sequence[-3].split("MEA=") [1]
	mea_line = mea_line.split("}")[0]
	final_line = my_sequence[-2]
	final_line = final_line.split("frequency of mfe structure in ensemble") [1]
	final_line = final_line.split(";")[0]
	ret_out = open(ret_url, "a")
	spam_writer = csv.writer(ret_out)
	genome_id = line_denuller(header.split("_"))
	ret = [genome_id[1:],header[1:], mea_line,final_line]
#	spam_writer.writerow(["Genome_id","Protein_id", "MEA", "Probability"])

	spam_writer.writerow(ret)
	return ret_list
