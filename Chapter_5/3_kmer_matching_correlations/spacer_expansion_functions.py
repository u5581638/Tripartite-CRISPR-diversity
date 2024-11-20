# spacer_expansion_functions
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

def kmer_pairwise_alignment(spacer,min_seq_id,coverage,spacer_num,spacer_id, contig, kmer_size=10):
	spacer_kmers = []
	i = 0
	spacer = spacer.upper()
	while (i + kmer_size < len(spacer)):
		new_kmer = spacer[i:i+kmer_size]
		my_kmer = Seq(new_kmer)
		the_kmer = SeqRecord(my_kmer,id=spacer_id + "_" +str(i))
		spacer_kmers.append(the_kmer)
		i += 1

	pairwise_alignments = []
	my_contig = str(contig.seq).upper()
	reverse_contig = reverse_complement(my_contig)
	switch = 0
	# may want to break and only allow one kmer match per spacer
	for kmer in spacer_kmers:
		m = re.search(str(kmer.seq),my_contig)
		mr = re.search(str(kmer.seq),reverse_contig)

		if (m):
			match_start = m.start()
			match_end = m.end()
			match_sequence = m.group(0)
			switch = 1
			break
		if (mr):
			match_start = len(contig) - mr.start()
			match_end = len(contig) - mr.end()
			match_sequence = reverse_complement(mr.group()) # take the reverse complement to be consistent with the normal frame
			switch = 1
			break
	if (switch == 1):
		outrow = [spacer_id, contig.id, len(match_sequence) / len(spacer),len(contig.seq),match_sequence,"NA","NA","NA",match_start,match_end,"NA","NA"]
		return outrow
	else:
		return None

# Do the pairwise kmer alignment as above, but record the size and positions of the kmer within the output matches.
def kmer_pairwise_alignment_query_length(spacer,min_seq_id,coverage,spacer_num,spacer_id, contig, kmer_size=10):
	spacer_kmers = []
	i = 0
	spacer = spacer.upper()
	while (i + kmer_size < len(spacer)):
		new_kmer = spacer[i:i+kmer_size]
		my_kmer = Seq(new_kmer)
		the_kmer = SeqRecord(my_kmer,id=spacer_id + "_" + "kstart:" + str(i) + "kend:" + str(i + kmer_size))
		spacer_kmers.append(the_kmer)
		i += 1

	pairwise_alignments = []
	my_contig = str(contig.seq).upper()
	reverse_contig = reverse_complement(my_contig)
	switch = 0
	# may want to break and only allow one kmer match per spacer
	for kmer in spacer_kmers:
		m = re.search(str(kmer.seq),my_contig)
		mr = re.search(str(kmer.seq),reverse_contig)
		query_len = len(spacer)
		qend = kmer.id.split("kend:") [0]
		qstart = qend.split("kstart:") [-1]
		qend = kmer.id.split("kend:") [-1]
		if (m):
			match_start = m.start()
			match_end = m.end()
			match_sequence = m.group(0)
			switch = 1
			break
		if (mr):
			# need to check the directions are correct
			match_start = len(contig) - mr.start()
			match_end = len(contig) - mr.end()
			match_sequence = reverse_complement(mr.group()) # take the reverse complement to be consistent with the normal frame
			switch = 1
			break
	if (switch == 1):
		outrow = [spacer_id, contig.id, len(match_sequence) / len(spacer),len(contig.seq),match_sequence,query_len,qstart,qend,match_start,match_end,"NA","NA"]
		return outrow
	else:
		return None

def kmer_pairwise_alignment_query_length_spacer_coord(spacer,min_seq_id,coverage,spacer_num,spacer_id, contig, kmer_size=10):
	spacer_kmers = []
	i = 0
	spacer = spacer.upper()
	spacer_size = len(spacer)
	while (i + kmer_size < len(spacer)):
		new_kmer = spacer[i:i+kmer_size]
		my_kmer = Seq(new_kmer)
		the_kmer = SeqRecord(my_kmer,id=spacer_id + "_" + "kstart:" + str(i) + "kend:" + str(i + kmer_size)) 
		spacer_kmers.append(the_kmer)
		i += 1

	pairwise_alignments = []
	my_contig = str(contig.seq).upper()
	reverse_contig = reverse_complement(my_contig)
	switch = 0
	# may want to break and only allow one kmer match per spacer
	for kmer in spacer_kmers:
		m = re.search(str(kmer.seq),my_contig)
		mr = re.search(str(kmer.seq),reverse_contig)
		query_len = len(spacer)
		qend = kmer.id.split("kend:") [0]
		qstart = qend.split("kstart:") [-1]
		qend = kmer.id.split("kend:") [-1]
		if (m):
			match_start = m.start() - int(qstart)
			match_end = m.end() + abs(spacer_size - int(qend) )
			match_sequence = m.group(0)
			switch = 1
			break
		if (mr):
			# need to check the directions are correct
			# Should be (contig) - mr.start() + abs(qstart)?
			match_start = len(contig) - mr.start() + abs(spacer_size - int(qstart)) # kmer start and end coordinates
			# Should be len(contig) - mr.end() - abs(spacer_size - qend)
			match_end = len(contig) - mr.end() - int(qend)
			# In the end the mipoint remains unchanged so the original math may be fine even though incorrect!!
			# May be useful to delete the reverse complement section so the matching sequence can be checked
			match_sequence = reverse_complement(mr.group(0)) # take the reverse complement to be consistent with the normal frame
			switch = 1
			break

	if (switch == 1):
		if (int(match_start) < 0 ):
			match_start = 0  
		if (int(match_end) < 0 ):
			match_end = 0 
		if (int(match_start) > len(contig)):
			match_start = len(contig)
		if (int(match_end) > len(contig)):
			match_end = len(contig)		
		outrow = [spacer_id, contig.id, len(match_sequence) / len(spacer),len(contig.seq),match_sequence,query_len,qstart,qend,match_start,match_end,"NA","NA"]
		return outrow
	else: 
		return None

def kmer_pairwise_alignment_regex(spacer,min_seq_id,coverage,spacer_num,spacer_id, contig, kmer_size=10):
	spacer_kmers = []
	i = 0
	spacer = spacer.upper()
	while (i + kmer_size < len(spacer)):
		new_kmer = spacer[i:i+kmer_size]
		my_kmer = Seq(new_kmer)
		the_kmer = SeqRecord(my_kmer,id=spacer_id + "_" +str(i))
		spacer_kmers.append(the_kmer)
		i += 1

	pairwise_alignments = []
	my_contig = str(contig.seq).upper()
	reverse_contig = reverse_complement(my_contig)
	switch = 0
	# may want to break and only allow one kmer match per spacer
	con_kmer = ""
	con_kmer = str(spacer_kmers[0].seq)
	i = 1
	while (i < len(spacer_kmers)):
		con_kmer += "|" + str(spacer_kmers[i].seq)
		i += 1
	con_pattern = re.compile(con_kmer)	
	m = re.search(con_pattern,my_contig)
	mr = re.search(con_pattern,reverse_contig)

	if (m):
		match_start = m.start()
		match_end = m.end()
		match_sequence = m.group(0)
		switch = 1
	if (mr):
		match_start = len(contig) - mr.start()
		match_end = len(contig) - mr.end()
		match_sequence = reverse_complement(mr.group()) # take the reverse complement to be consistent with the normal frame
		switch = 1
	if (switch == 1):
		outrow = [spacer_id, contig.id, len(match_sequence) / len(spacer),len(contig.seq),match_sequence,"NA","NA","NA",match_start,match_end,"NA","NA"]
		return outrow
	else:
		return None

def blast_pairwise_alignment(spacer,min_seq_id,coverage,spacer_num,spacer_id):
	# write spacer to file
	my_spacer = Seq(spacer)
	the_spacer = SeqRecord(my_spacer,id=spacer_id)
	SeqIO.write(the_spacer,"spacer1.fasta","fasta")
	
	subprocess.run(["blastn -db " + "contig1.fasta" + " -query " + "spacer1.fasta" + " -outfmt 10 " + " -out " + "pairwise_aligned.csv" +  " -max_target_seqs 1000000" + " -evalue 0.05 -word_size 4 -qcov_hsp_perc 50 -perc_identity 50" ],shell=True)
	try:
		with open("pairwise_aligned.csv") as csvfile:
			pairwise_alignments = list(csv.reader(csvfile))
	except:
		print("Caution:")
		return None	
	# may want to comment out this line in testing
	subprocess.run("rm spacer1.fasta pairwise_aligned.csv",shell=True)
	ret_row = []
	for row in pairwise_alignments:
		ret_row.append([row[2],row[8],row[9],spacer])
	print("End of BLAST:")
#	print(pairwise_alignments)
	if (pairwise_alignments == []):
		return None
	elif (pairwise_alignments[0][1] == "genome_block_19212|Ga0315285_10012069:1-11070" and pairwise_alignments[0][8] == "958" and pairwise_alignments[0][9] == "989"):
		print("Done!!")
	else:	
		return pairwise_alignments

# replace target site of original mapped spacer in contig by "N"
def id_sorter(x):
	return x[-1]
def spacer_sorter(x):
	return x[1]

def scoring(seq1,seq2):
# sequences must be the same length for this to work
	i = 0
	my_score = 0
	while (i < len(seq1)):
		if (seq1[i] == '-' and seq2[i] == '-'):
			i += 1
			continue
		elif(seq1[i] == seq2[i]):
			my_score += 1
			i += 1
		else:
			i += 1
	return my_score

def mask_contig (phage_contig, start, end):
	phage_seq = str(phage_contig.seq)
	if (start < end):
		phage_seq_start = phage_seq[:start]
		phage_seq_end = phage_seq[end:]
		phage_diff = end - start
	else:
		phage_seq_start = phage_seq[:end]
		phage_seq_end = phage_seq[start:]
		phage_diff = start - end
	insert_seq = "N" * phage_diff
	phage_seq = Seq(phage_seq_start + insert_seq + phage_seq_end)
	# I hope this assignment works for parameters like a normal variable
	phage_contig.seq = phage_seq
	return phage_contig

def coords_in_range(coords, phage_start,phage_end):
	ret_coords = []
	for coord in coords:
		if (phage_start <= coord[0] <= phage_end and phage_start <= coord[1] <= phage_end):
			ret_coords.append(coord)
	return ret_coords


# eliminate spacers with homology to original spacer + original spacer.
def identity(spacer1,spacer2):
	spacer1_id = spacer1.id 
	spacer2_id = spacer2.id
	alignment = pairwise2.align.globalxx(spacer1.seq,spacer2.seq)
	spacer1.seq = alignment[0][0]  
	spacer2.seq = alignment[0][1]
	i = 0
	count = 0
	return [spacer1_id,spacer2_id,alignment[0][2]/len(spacer1.seq)]

def spacer_identity(new_arr, the_spacers):
	ret_arrs = []
	for arr in new_arr:
		for spacer in the_spacers:
			hit = identity(spacer,arr) # [arr_id,spacer_id,perc_id]
			if (hit[2] > 0.3):
				ret_arrs.append(hit)
	return ret_arrs

def gen_synthetic_phage(phage_length,phage_id):

	ret_str = []
	i = 0
	nucleotides = ["A","T","C","G"]
	while (i < phage_length):
		my_nucleotide = random.choice(nucleotides)
		ret_str.append(my_nucleotide)
		i += 1
	return SeqRecord(Seq("".join(ret_str)),id=phage_id)



def kmer_homolog_elimination(array,mapped_spacers):
	arr_spacers = []
	new_arr = []
	for row in array:
		spacer_id = row[0] + "|" + "spacer_start_pos:" + str(row[8]) + "|" + "spacer_end_pos:" + str(row[9]) + "|" + "global_start_pos:" + str(int(row[8])) + "|" + "global_end_pos:" + str(int(row[9])) + "|" + "array_tool:" + row[15]
		the_spacer = SeqRecord(Seq(row[14]),id=spacer_id)
		arr_spacers.append(the_spacer)
		new_arr.append([spacer_id] + row[1:])
	the_spacers = []	
	for row in mapped_spacers:
		the_spacer = SeqRecord(Seq(row[26]),id=row[0])
		the_spacers.append(the_spacer)

	hit_table = spacer_identity(arr_spacers,the_spacers)
	mapping_ids = set()
	for hit in hit_table:
		if (hit[1] not in mapping_ids):
			mapping_ids.add(hit[1])

	ret_arrs = []
	for arr in new_arr:
		if (arr[0] not in mapping_ids):
			ret_arrs.append(arr)

	return ret_arrs




def blast_homolog_elimination(array,mapped_spacers):
	# BLAST mapped spacers and eliminate the mapped spacers and any homologs
	arr_spacers = []
	new_arr = []
	for row in array:
		spacer_id = row[0] + "|" + "spacer_start_pos:" + str(row[8]) + "|" + "spacer_end_pos:" + str(row[9]) + "|" + "global_start_pos:" + str(int(row[8])) + "|" + "global_end_pos:" + str(int(row[9])) + "|" + "array_tool:" + row[15]
		the_spacer = SeqRecord(Seq(row[14]),id=spacer_id)
		arr_spacers.append(the_spacer)
		new_arr.append([spacer_id] + row[1:])
		
	SeqIO.write(arr_spacers, "spacer_db.fasta","fasta")
	subprocess.run(["makeblastdb -in " + "spacer_db.fasta" + " -dbtype nucl"],shell=True)
	the_spacers = []
	for row in mapped_spacers:
		the_spacer = SeqRecord(Seq(row[26]),id=row[0])
		the_spacers.append(the_spacer)
	SeqIO.write(the_spacers,"mapped_spacers.fasta","fasta")
	subprocess.run(["blastn -db " + "spacer_db.fasta" + " -query " + "mapped_spacers.fasta" + " -outfmt 10 " + " -out " + "mapped_spacers.csv" + " -perc_identity 50 " +  " -max_target_seqs 1000000" + " -evalue 0.05 -word_size 4 -qcov_hsp_perc 50 " ],shell=True)
	with open("mapped_spacers.csv","r") as csvfile:
		hit_table = list(csv.reader(csvfile))
		mapping_ids = set()
		for hit in hit_table:
			if (hit[1] not in mapping_ids):
				mapping_ids.add(hit[1])
		ret_arrs = []
		

		for arr in new_arr:
		#	print(arr[0])
			if arr[0] not in mapping_ids:
				ret_arrs.append(arr)
				print("yay!!")
			else:
				print("Self_target_check!")

	subprocess.run(["rm mapped_spacers.fasta spacer_db.fasta mapped_spacers.csv"],shell=True)
	return ret_arrs
