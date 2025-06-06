# mapped_spacer_expansion

# map spacers onto contigs for which a pre-existing match already exists
# This requires a multi-step approach.

# 1. Need to develop a pairwise alignment tool to substitute BLAST, which is unsuitable due to it's need to be indexed as a database.

# 2. Need to map all spacers excluding the original mapped spacer used to identify the contig. Should also mask the original mapped region as 
# 	 There may be duplicates in the CRISPR array encoding the same spacer.
# 3. Need to compare the result of mapping at different thresholds with the hit rate from mapping against a randomised db of similar number and size of sequences (in bytes)
# 4. It may be best if these randomised sequences are sourced from the same samples/genome blocks to show the effect is not sample specific.

# 1. development of a pairwise tool to determine an identity score, identify all S-W matches, and exclude matching targets below this threshold. 
#    The matched target start and end coordinates, as well as other properties should be saved and used to build a csv file in the same format as BLAST output 10.
#	 This will simplify subsequent code to find primed spacers.

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
import spacer_expansion_functions

# INPUT: 1. table containing the concensous arrays predicted and validated (see chapter 4) for a given subtype. 
#		 i.e. cas12a.fasta_all_hits.csv_genomes.fasta_crisprs.lst_reconciled_full_arr_positions.csv_no_blanks.csv
#		 2. spacer hitmap table giving the spacer hits to target sequences
#		 i.e. cas12a.fasta_all_hits.csv_genomes.fasta_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_hitmap.csv_standardised.csv_non_redundant.csv_filtered.csv
#		 3. a file containing the target sequence (phage) contigs
#		 i.e. cas12a.fasta_all_hits.csv_genomes.fasta_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_hitmap.csv_standardised.csv_non_redundant.csv_filtered.csv_deduplicated.csv_filtered_hits_extracted_faidx_bp_window.fasta
# OUTPUT: 1. table containing kmer matches on each contig
#		  i.e cas12a_detection_parallel_gadi_kmer14_29-3-2024.csv
#		  2. table containing kmer matches to associated direct repeats on each contig
#		  i.e. cas12a_detection_parallel_gadi_kmer14_29-3-2024.csv_dr.csv
# SHELL: see cas12a_kmer14_detection.sh

kmer_switch = 1
with open(sys.argv[1], "r") as csvfile:
	arrays = list(csv.reader(csvfile))

# take largest concensus array and number spacers from PPS
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

# next need to load mapped spacer genes
with open(sys.argv[2],"r") as csvfile:
	dedup_spacer_hits = list(csv.reader(csvfile))

dedup_spacer_dict = {}
phage_id_set = {}
# initialise a dictionary of mapped phage hits by crispr array as well as a non-redundant set of phage_ids
for spacer_hit in dedup_spacer_hits[1:]:
	if (spacer_hit[1] not in phage_id_set):
		phage_id_set[spacer_hit[1]] = [[spacer_hit[8],spacer_hit[9]]]
	else:
		phage_id_set[spacer_hit[1]].append([spacer_hit[8],spacer_hit[9]])	
	genome_id = spacer_hit[0].split("|")[0]
	array_no = spacer_hit[-2]
	# double check this key works correctly
	genome_id_arr = (genome_id, str(array_no))
	if (genome_id_arr) not in dedup_spacer_dict:
		dedup_spacer_dict[genome_id_arr] = [spacer_hit]
	else:
		dedup_spacer_dict[genome_id_arr].append(spacer_hit)
dedup_spacer_list = list(dedup_spacer_dict.values()) # may be best iterating on the dictionary directly

mapped_phages = list(SeqIO.parse(sys.argv[3], "fasta"))
mapped_phage_dict = {}

# Impossible to map phages directly due to being segments. Instead made in to a dict mapping to a list of phage contigs from the same genome. Will then need to iterate over the range of phages
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
# will actually need to mask every mapped sequence in the mapped phage for effectiveness. This will require dictionalising all mapped hits which target the same phage sequence
i = 0
phage_count = 0

# masking all mapped_phage_contigs:
# 1. iterate through the mapped spacer list
# 2. Group together all the same phage ids
# 3. Retrieve the phage contig windows
# 4. Mask all existing target sites.
# 5. Assign the new masked contigs to the new mapped_phage_dict.

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

# create a set of all mapped_phage_ids
# need to mask all mapped contigs before doing the main loop:
# start main loop:
# In this loop:
# 1. Iterate through mapped sequences.
# 2. Retrieve mapped phage sequences with coordinates overlapping the target range. A single contig should be sufficent!!
# 3. Deduplicate homologous spacers
# 4. Perform pairwise alignments
# 5. Save the original mapped alignment + additional alignments in the same mapped format as prior used to PPS calculation (may have to put NA for some columns).
# 6. Write these entries to a file in append/write mode.
# mapped phages need to be dictionalised so that the appropiate corresponding mapped gene can be looked up

for array in dedup_spacer_list:
	# group spacers in array into lists based on a conserved phage contig`
	phage_id_dict = {}
	for spacer in array:
		i += 1
		if (spacer[1] not in phage_id_dict):
			phage_id_dict[spacer[1]] = [spacer]
		else:
			phage_id_dict[spacer[1]].append(spacer)
	
	matching_arr = array_dict[(spacer[0].split("|")[0],int(spacer[29]))]
	for phage_id in phage_id_dict:
		# need to sort according to mapped phage genome.
		ms_rows = []
		for row in phage_id_dict[phage_id]:
			ms_rows.append(row) # Is there a dictionalised alternative to using ms_rows
			phage_count += 1	
			 # In theory this should never return a key error, as all spacers come from this file	
				# code to eliminate spacers with homology to the existing spacers using BLAST and target coord info
		if (kmer_switch != 1):
			matching_arr = spacer_expansion_functions.blast_homolog_elimination(matching_arr, ms_rows)
		else:
			matching_arr = spacer_expansion_functions.kmer_homolog_elimination(matching_arr,ms_rows)

		phage_contigs = mapped_phage_dict[phage_id] # This should never return an error if the input files are correct
		# end of code to eliminate homologous spacers
		for new_phage_contig in phage_contigs:
			phage_contig_start = new_phage_contig.id.split(":") [1]
			phage_contig_end = phage_contig_start.split("-") [1]
			phage_contig_start = phage_contig_start.split("-") [0]



			original_spacers = []
			if (kmer_switch != 1):
				SeqIO.write(new_phage_contig, "contig1.fasta","fasta")
				subprocess.run(["makeblastdb -in " + "contig1.fasta" + " -dbtype nucl"],shell=True)
			# Exclude the original spacers + any spacers with homology
			for arr_spacer in matching_arr:
				# should be able to eliminate this as these targets should be masked.
				if (kmer_switch == 1):
					# for each spacer in the array (bar the original) do the kmer search
					pairwise_mapping = spacer_expansion_functions.kmer_pairwise_alignment_query_length_spacer_coord(arr_spacer[14],20,5,arr_spacer[-1],arr_spacer[0], new_phage_contig,10)
					# kmer search using Direct repeats with spacer kmer handles.
					pairwise_mapping2 = spacer_expansion_functions.kmer_pairwise_alignment_query_length_spacer_coord(arr_spacer[14][-5:] + arr_spacer[12] + arr_spacer[14][:5],20,5,arr_spacer[-1],arr_spacer[0], new_phage_contig,10)
				else:
					pairwise_mapping = spacer_expansion_functions.blast_pairwise_alignment(arr_spacer[14], 20,5, arr_spacer[-1],arr_spacer[0])
				if (pairwise_mapping is not None):
					if (kmer_switch != 1):
						pairwise_mapping = pairwise_mapping[0]
					outrow = []
				# convert to the arr_spacer_equivalent. Need to consider the reverse complement!!
					target_start = str( int(pairwise_mapping [8])  + int(phage_contig_start) - 1)
					target_end = str( int(pairwise_mapping [9]) + int(phage_contig_start) - 1)
					start_index = arr_spacer[0].split("::")[1]
					start_index = start_index.split(":")[0]				
					mapped_spacer_id = arr_spacer[0].split("|")[0] + "|" + "spacer_start_pos:" + str(arr_spacer[8]) + "|" + "spacer_end_pos:" + str(arr_spacer[9]) + "|" + "global_start_pos:" + str(int(arr_spacer[8]) + int(start_index) - 1) + "|" + "global_end_pos:" + str(int(arr_spacer[9]) + int(start_index) - 1) + "|" + "array_tool:" + arr_spacer[15]
					spamwriter.writerow([mapped_spacer_id] + [pairwise_mapping[1].split(":")[0]] + pairwise_mapping[2:8] + [target_start,target_end] + ["NA","NA","NA"] + [arr_spacer[1]] + ["NA","NA","NA","NA", arr_spacer[6],arr_spacer[7],arr_spacer[8],arr_spacer[9],arr_spacer[10],arr_spacer[11],arr_spacer[12],arr_spacer[13],arr_spacer[14],arr_spacer[15],arr_spacer[16],arr_spacer[-2],arr_spacer[-1]])
				if (pairwise_mapping2 is not None):
					if (kmer_switch != 1):
						pairwise_mapping2 = pairwise_mapping2[0]
					outrow = []
				# convert to the arr_spacer_equivalent. Need to consider the reverse complement!!
					target_start = str( int(pairwise_mapping2 [8])  + int(phage_contig_start) - 1)
					target_end = str( int(pairwise_mapping2 [9]) + int(phage_contig_start) - 1)
					start_index = arr_spacer[0].split("::")[1]
					start_index = start_index.split(":")[0]				
					mapped_spacer_id = arr_spacer[0].split("|")[0] + "|" + "spacer_start_pos:" + str(arr_spacer[8]) + "|" + "spacer_end_pos:" + str(arr_spacer[9]) + "|" + "global_start_pos:" + str(int(arr_spacer[8]) + int(start_index) - 1) + "|" + "global_end_pos:" + str(int(arr_spacer[9]) + int(start_index) - 1) + "|" + "array_tool:" + arr_spacer[15]
					spamwriter2.writerow([mapped_spacer_id] + [pairwise_mapping2[1].split(":")[0]] + pairwise_mapping2[2:8] + [target_start,target_end] + ["NA","NA","NA"] + [arr_spacer[1]] + ["NA","NA","NA","NA", arr_spacer[6],arr_spacer[7],arr_spacer[8],arr_spacer[9],arr_spacer[10],arr_spacer[11],arr_spacer[12],arr_spacer[13],arr_spacer[14],arr_spacer[15],arr_spacer[16],arr_spacer[-2],arr_spacer[-1]])
			if (kmer_switch != 1):
				subprocess.run("rm contig1.fasta",shell=True)
					
print("phage_count:")
print(phage_count)	
ret_out.close()
ret_out2.close()

# With control sequence the genome to genome comparison should be random - no mapped genome as reference. Should be against a block of roughly the same size in bp