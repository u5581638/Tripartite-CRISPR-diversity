# mapped_spacer_expansion_randomised_phage

# script to perform a control local pairwise alignment by reshuffling the mapped phages.

# main difference will be that the mapped phages are only used to retrieve the set of possible arrays.
# each array is then assigned at random without duplication or replacement to each phage contig.
# This means that filtered by spacer homology is likely not nessessary in theory because spacers should not map to the non-primary array.
# May be wise to filter the matching spacer anyway though to prevent hits to homologous contigs if there are multiple copies of the same phage (which could artificallu inflate the results.)
# masking a contig also may not make sense.

from Bio import Align
from Bio import SeqIO
import csv
import sys
from Bio.Seq import reverse_complement
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import filtered_arrs_appender_SQL4
import re
import copy
import subprocess
import random
import spacer_expansion_functions

# split spacer into kmers and match these kmers against the target genome.
# This should be significantly faster than BLAST pairwise alignments
# need to format the output to match the existing ()

kmer_switch = 1
with open(sys.argv[1], "r") as csvfile:
	arrays = list(csv.reader(csvfile))

arrays = filtered_arrs_appender_SQL4.largest_array_merge(arrays)
arrays = filtered_arrs_appender_SQL4.number_arrays(arrays[1:])
arr_out = open(sys.argv[1] + "arr_cpy.csv","w")
spamwriter = csv.writer(arr_out)
for arr in arrays:
	for row in arr:
		spamwriter.writerow(row)
arr_out.close()

array_dict = {}
for arr_group in arrays:
	for arr in arr_group:
		arr_id = (arr[0],arr[-2]) # array number should be second from the right
		if (arr_id not in array_dict):
			array_dict[arr_id] = [arr]
		else:
			array_dict[arr_id].append(arr)	
# code for array dictionary goes below:

# next need to load mapped spacer genes and mask
with open(sys.argv[2],"r") as csvfile:
	dedup_spacer_hits = list(csv.reader(csvfile))

dedup_spacer_dict = {}
phage_id_set = {}
i = 0
phage_count = 0
for spacer_hit in dedup_spacer_hits[1:]:
	if (spacer_hit[1] not in phage_id_set):
		phage_id_set[spacer_hit[1]] = [[spacer_hit[8],spacer_hit[9]]]
	else:
		phage_id_set[spacer_hit[1]].append([spacer_hit[8],spacer_hit[9]])

	genome_id = spacer_hit[0].split("|")[0]
	array_no = spacer_hit[-2]
	genome_id_arr = (genome_id, str(array_no))
	if (genome_id_arr) not in dedup_spacer_dict:
		dedup_spacer_dict[genome_id_arr] = [spacer_hit]
	else:
		dedup_spacer_dict[genome_id_arr].append(spacer_hit)
dedup_spacer_list = list(dedup_spacer_dict.values()) # may be best iterating on the dictionary directly

mapped_phages = list(SeqIO.parse(sys.argv[3], "fasta"))
mapped_phage_dict = {}

for m_phage in mapped_phages:
	phage_id = m_phage.id.split(":") [0]
	if (phage_id not in mapped_phage_dict):
		mapped_phage_dict[phage_id] = [m_phage]
	else:
		mapped_phage_dict[phage_id].append(m_phage)

ret_out = open(sys.argv[4],"a")
spamwriter = csv.writer(ret_out)
# write the header
spamwriter.writerow(["Spacer_id","Phage_id","Perc_id","Length", "Mismatches","Gapopen","query_start","query_end","Mapped_start_site","Mapped_end_site","evalue","bitscore","Genome_id","orientation","orientation_score","orientation_confidence","questionable_array","array_score","CRISPR-start","CRISPR-end","repeat_start","repeat_end","spacer_start","spacer_end","dr_repeat_original","dr_repeat_concensous","spacer","Array_tool","array_number","spacer_number"])

ret_out2 = open(sys.argv[4] + "_dr.csv","a")
spamwriter2 = csv.writer(ret_out2)
# write the header
spamwriter2.writerow(["Spacer_id","Phage_id","Perc_id","Length", "Mismatches","Gapopen","query_start","query_end","Mapped_start_site","Mapped_end_site","evalue","bitscore","Genome_id","orientation","orientation_score","orientation_confidence","questionable_array","array_score","CRISPR-start","CRISPR-end","repeat_start","repeat_end","spacer_start","spacer_end","dr_repeat_original","dr_repeat_concensous","spacer","Array_tool","array_number","spacer_number"])
# want to select just the subset of arrays corresponding to the mapped spacers.
for phage_set_id in phage_id_set:
	i = 0
	while (i < len(mapped_phage_dict[phage_set_id])):
		phage_contig_start = mapped_phage_dict[phage_set_id][i].id.split(":") [1]
		phage_contig_end = phage_contig_start.split("-") [1]
		phage_contig_start = phage_contig_start.split("-") [0]
		target_coords = spacer_expansion_functions.coords_in_range(phage_id_set[phage_set_id],phage_contig_start,phage_contig_end)
		if (len(target_coords) > 0):
			for coord in target_coords:
				mapped_phage_dict[phage_set_id][i] = spacer_expansion_functions.mask_contig(mapped_phage_dict[phage_set_id][i], int(coord[0]) - int(phage_contig_start) + 1,int(coord[1]) - int(phage_contig_start) + 1)
		i += 1

phage_pool = []

for array in dedup_spacer_list:
	phage_id_dict = {}
	for spacer in array:
		print("i:")
		print(i)
		i += 1
		if (spacer[1] not in phage_id_dict):
			phage_id_dict[spacer[1]] = [spacer]
		else:
			phage_id_dict[spacer[1]].append(spacer)
	for phage_id in phage_id_dict:
		phage_pool.extend(mapped_phage_dict[phage_id])


for array in dedup_spacer_list:
	# group spacers in array into lists based on a conserved phage contig`
	phage_id_dict = {}
	for spacer in array:
		print("i:")
		print(i)
		i += 1
		if (spacer[1] not in phage_id_dict):
			phage_id_dict[spacer[1]] = [spacer]
		else:
			phage_id_dict[spacer[1]].append(spacer)
	# Iterate over the phage ids from the mapped spacers
	matching_arr = array_dict[(spacer[0].split("|")[0],int(spacer[29]))]
	for phage_id in phage_id_dict:
		ms_rows = []
		for row in phage_id_dict[phage_id]:
			ms_rows.append(row)
		
		if (kmer_switch != 1):
			matching_arr =  spacer_expansion_functions.blast_homolog_elimination(matching_arr, ms_rows)
		else:
			matching_arr = spacer_expansion_functions.kmer_homolog_elimination(matching_arr,ms_rows)
		# take a random phage key
		contig_number = len(mapped_phage_dict[phage_id])
		mapped_phage_random = []
		k = 0 
		while (k < contig_number):
			phage_index = random.randint(0,len(phage_pool) - 1)
			phage_contigv = phage_pool.pop(phage_index)
			mapped_phage_random.append(phage_contigv)
			k += 1
		# take 4 random contigs from the pool of phage-mapped contigs without replacement.

		 # should this be without replacement??
		for new_phage_contig in mapped_phage_random:	
			phage_contig_start = new_phage_contig.id.split(":") [1]
			phage_contig_end = phage_contig_start.split("-") [1]
			phage_contig_start = phage_contig_start.split("-") [0]	
			# this is the first line that needs to change.
			original_spacers = []
			phage_count += 1
			if (kmer_switch != 1):
				SeqIO.write(new_phage_contig, "contig1.fasta","fasta")
				subprocess.run(["makeblastdb -in " + "contig1.fasta" + " -dbtype nucl"],shell=True)
			# Exclude the original spacers + any spacers with homology (use BLAST) -> could use .
			for arr_spacer in matching_arr:
					# should be able to eliminate this as these targets should be masked.
					# should include the option to do kmer comparison instead.
				if (kmer_switch == 1):
					pairwise_mapping = spacer_expansion_functions.kmer_pairwise_alignment_query_length_spacer_coord(arr_spacer[14],20,5,arr_spacer[-1],arr_spacer[0], new_phage_contig,10)
					pairwise_mapping2 = spacer_expansion_functions.kmer_pairwise_alignment_query_length_spacer_coord(arr_spacer[14][-5:] + arr_spacer[12] + arr_spacer[14][:5],20,5,arr_spacer[-1],arr_spacer[0], new_phage_contig,10)

				else:
					pairwise_mapping = spacer_expansion_functions.blast_pairwise_alignment(arr_spacer[14], 20,5, arr_spacer[-1],arr_spacer[0],target_coords)
				if (pairwise_mapping is not None):
					if (kmer_switch != 1):
						pairwise_mapping = pairwise_mapping[0]
					outrow = []
					target_start = str( int(pairwise_mapping [8])  + int(phage_contig_start) - 1)
					target_end = str( int(pairwise_mapping [9]) + int(phage_contig_start) - 1)
					start_index = arr_spacer[0].split("::")[1]
					start_index = start_index.split(":")[0]				
					mapped_spacer_id = arr_spacer[0].split("|")[0] + "|" + "spacer_start_pos:" + str(arr_spacer[8]) + "|" + "spacer_end_pos:" + str(arr_spacer[9]) + "|" + "global_start_pos:" + str(int(arr_spacer[8]) + int(start_index) - 1) + "|" + "global_end_pos:" + str(int(arr_spacer[9]) + int(start_index) - 1) + "| " + "array_tool:" + arr_spacer[15]
					# add evalue and scores
					spamwriter.writerow([mapped_spacer_id] + [pairwise_mapping[1].split(":")[0]] + pairwise_mapping[2:8] + [target_start,target_end] + pairwise_mapping[10:11] + ["NA","NA","NA","NA","NA","NA", arr_spacer[6],arr_spacer[7],arr_spacer[8],arr_spacer[9],arr_spacer[10],arr_spacer[11],arr_spacer[12],arr_spacer[13],arr_spacer[14],arr_spacer[15],arr_spacer[16],arr_spacer[-2],arr_spacer[-1]])
				if (pairwise_mapping2 is not None):
					if (kmer_switch != 1):
						pairwise_mapping2 = pairwise_mapping2[0]
					outrow = []
					target_start = str( int(pairwise_mapping2 [8])  + int(phage_contig_start) - 1)
					target_end = str( int(pairwise_mapping2 [9]) + int(phage_contig_start) - 1)
					start_index = arr_spacer[0].split("::")[1]
					start_index = start_index.split(":")[0]				
					mapped_spacer_id = arr_spacer[0].split("|")[0] + "|" + "spacer_start_pos:" + str(arr_spacer[8]) + "|" + "spacer_end_pos:" + str(arr_spacer[9]) + "|" + "global_start_pos:" + str(int(arr_spacer[8]) + int(start_index) - 1) + "|" + "global_end_pos:" + str(int(arr_spacer[9]) + int(start_index) - 1) + "| " + "array_tool:" + arr_spacer[15]
					# add evalue and scores
					spamwriter2.writerow([mapped_spacer_id] + [pairwise_mapping2[1].split(":")[0]] + pairwise_mapping2[2:8] + [target_start,target_end] + pairwise_mapping2[10:11] + ["NA","NA","NA","NA","NA","NA", arr_spacer[6],arr_spacer[7],arr_spacer[8],arr_spacer[9],arr_spacer[10],arr_spacer[11],arr_spacer[12],arr_spacer[13],arr_spacer[14],arr_spacer[15],arr_spacer[16],arr_spacer[-2],arr_spacer[-1]])
			if (kmer_switch != 1):
				subprocess.run("rm contig1.fasta",shell=True)
							# consider some control statements to flag errors in case things don't.
					# for each spacer in the array (bar the original) do the pairwise alignment
				# The above if statement should always be triggered at least once unless something is broken with the mapped genome
	# can I map arrays to their mapped phages by id? Use tuples as key?

print("phage_count:")
print(phage_count)
phage_count += 1	
ret_out.close()
ret_out2.close()



