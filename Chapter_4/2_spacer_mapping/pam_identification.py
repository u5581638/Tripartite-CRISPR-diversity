# program to find PAM sequences from mapped genomes and genome_hit files and write upstream and downstream web logos as well as a tabulated copy of the upstream and downstream potential PAM sequences for each mapped hit
from Bio import SeqIO
import sys
import csv
from Bio import AlignIO
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord

def rev_comp(genome):
	ret_str = ""
	for ele in genome:
		if (ele == "A" or ele == "a"):
			ret_str += "T"
		elif (ele == "T" or ele == "t"):
			ret_str += "A"
		elif (ele == "C" or ele == "c"):
			ret_str += "G"
		elif (ele == "G" or ele == "g"):
			ret_str += "C"
		else:
			ret_str += ele
	ret_str = ret_str[::-1]		
	return ret_str						

def slice (mapped_genome, start, end, offset):
#	print("slice values:")
#	print(start,end, len(mapped_genome))
	if (start - offset < 0 and end + offset < len(mapped_genome)):
		return (mapped_genome[:start], mapped_genome[end:(end + offset)])
	elif (start - offset >= 0 and end + offset < len(mapped_genome)):
		return (mapped_genome[start - offset:start], mapped_genome[end:(offset+end)])
	elif (start - offset >= 0 and end + offset >= len(mapped_genome)):
		return (mapped_genome[start - offset:start], mapped_genome[end:])
	else:
		
		print("Warning, this should not happen, spacers too small!!!!!!!")
		return (mapped_genome[:start], mapped_genome[end:])	


 # Upstream and downstream PAMs must be aligned based on the CRISPR-orientation.
def pam_extraction(mapped_genome, spacer_start, spacer_end, blast_direction, crispr_direction, offset): # funciton
	offset = int(offset)
	if ((crispr_direction == "Forward" and blast_direction == 1) or (crispr_direction == "Reverse" and blast_direction == -1)): # get gaetan to check this logic!!
		return slice(str(mapped_genome.seq), spacer_start, spacer_end, offset)
	else:
		return slice(rev_comp(mapped_genome), len(mapped_genome) - spacer_end +1,len(mapped_genome) - spacer_start +1, offset)	
# 	return (upstream_pam, downstream_pam)

# 11
def pam_id (master_table_url, genomes_url, out_dir, offset,index):
	master_hit_handle = open(master_table_url, "r")
	master_hit = csv.reader(master_hit_handle)
	mapped_genomes = SeqIO.parse(genomes_url, "fasta")
	genome_dict = {}
	ret_out = open (out_dir, "w")
	spam_writer = csv.writer(ret_out)
	# need to modify this so that the spacer id is included!!!
	spam_writer.writerow(["Spacer_id", "Genome_id", "Phage_id", "Spacer_start", "spacer_end","crispr_direction","blast_direction","mapped_spacer_start","mapped_spacer_end","PAM_start","PAM_end","RUN"])
	full_out = open (master_table_url + "_w_pams.csv", "w")
	full_writer = csv.writer(full_out)
	print("Initialise!")
	genome_slice_coord_dict = {}
	for genome in mapped_genomes:
		genome_id = genome.id 
		genome_id = genome_id.split(":")
		genome_start = int(genome_id[1].split("-")[0]) - 1
		genome_end = int(genome_id[1].split("-")[1])
		genome_id = genome_id[0]
		# Haven't considered the case for multiple windows from the same contig id!!
		if (genome_id not in genome_dict):
			genome_dict[genome_id] = [[genome, genome_start, genome_end]]
		else:
			genome_dict[genome_id].append([genome, genome_start, genome_end])	
	#	genome_slice_coord_dict[genome_id] = genome_start # This might be important for figuring out where a given slice is meant to be!!

	master_hit_arrays = {}
	next(master_hit)
	header = ["Spacer_id","Phage_id","Perc_id","Length","Mismatches","Gapopen","query_start","query_end","Mapped_start_site","Mapped_end_site","evalue","bitscore","Genome_id","orientation","orientation_score","orientation_confidence","questionable_array","array_score","CRISPR-start","CRISPR-end","repeat_start","repeat_end","spacer_start","spacer_end","dr_repeat_original","dr_repeat_concensous","spacer","Array_tool","run","array_number","spacer_number"] + ["PAM_start", "PAM_end", "RUN"]
	full_writer.writerow(header)
	print("dict_matching!!")
	for hit in master_hit:
		my_id = hit[0] 
		if (my_id  not in master_hit_arrays):
			master_hit_arrays[my_id] = [hit]
		else:
			master_hit_arrays[my_id].append(hit)
	ret_PAM = []	
#	print(master_hit_arrays.values())
	print("dict_initialised!!")
	for array in master_hit_arrays.values():
 	
		for hit in array:
			if (int(hit[8]) > int(hit[9])):
				spacer_start = int(hit[9])
				spacer_end = int(hit[8])
				blast_direction = -1
			else:
				spacer_start = int(hit[8])
				spacer_end =  int(hit[9])
				blast_direction = 1
			
			crispr_direction = str(hit[13])
			# This has to be the critical line!
			# need a more precise means for selecting the start coord for each contig!!
			for entry in genome_dict[hit[1]]:
				if (entry[1] <= spacer_start <= entry[2] or entry[1] <= spacer_end <= entry[2]):
					genome_slice_coord = entry
					break

			pams = pam_extraction(genome_slice_coord[0],spacer_start - genome_slice_coord[1], spacer_end - genome_slice_coord[1], blast_direction, crispr_direction, offset)

			if (pams[0] == "" or pams[1] == ""):
				print("special_PAM")

				print(hit)
			ret_PAM.append(pams)
			id_start = hit[22] # ? spacer_start_coord_in_host genome
			id_end =  hit[23] # These are placeholders. Replace with the actual genome start ids
			out_row = [hit[0].split(" ")[0],hit[0].split("|")[0], hit[1], id_start, id_end, str(crispr_direction), str(blast_direction), hit[8], hit[9], str(pams[0]), str(pams[1]), index]
			full_row = hit + [pams[0]] + [pams[1]] + [index]
			spam_writer.writerow(out_row)
			full_writer.writerow(full_row)
 			# need to consider the best output format!!


	ret_out.close()
	full_out.close()
	master_hit_handle.close()
	return 0

