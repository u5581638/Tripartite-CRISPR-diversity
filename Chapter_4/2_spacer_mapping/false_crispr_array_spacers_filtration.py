# program to filter spacers which are actually entire crispr_arrays generated from PILER-CR/CRT
# group the spacers by id then delete any spacer more than 1.5 times the average size.


from Bio import SeqIO
import statistics

def filtration(input_url):
	spacers = SeqIO.parse(input_url, "fasta")
	seq_dict = {}

	ret_list = []
	for spacer in spacers:
		my_id = spacer.id.split("|") [0] # see if this works, it may not.
		if my_id not in seq_dict:
			seq_dict[my_id] = [spacer]
		else:
			seq_dict[my_id].append(spacer)

	for spacers in seq_dict.values():
		seq_lens = []
		for spacer in spacers:
			seq_lens.append(len(spacer.seq))
		avg_len = float(statistics.mean(seq_lens))
		length_limit_upper = 3 * avg_len
		length_limit_lower = 0.5 * avg_len
		for spacer in spacers: # need to check I can run this loop again. # Is the pointer in the right spot (at the start)?
			if (len(spacer.seq) < length_limit_upper and len(spacer.seq) > length_limit_lower):
				ret_list.append(spacer)
	SeqIO.write(ret_list, input_url + "_arrs_filtered.fasta", "fasta")
	return 0


