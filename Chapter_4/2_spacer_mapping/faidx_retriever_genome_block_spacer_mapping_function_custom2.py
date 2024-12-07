# program to utilise samtools faidx to retrieve genomes
import csv
import sys
import copy
import subprocess
from Bio import SeqIO

# In this case the program should take an input genome/csv entry and load the genome based on the block id, then use the index (genome name) to retrieve the exact genome.

# script to retrieve mapped sequence genomes by identifiers using samtools.

def mapped_genome_retrieval(spacer_mapping_table_url, extraction_range=20000):
	csvfile = open(spacer_mapping_table_url, "r")
	hit_table = csv.reader(csvfile)

	extraction_range = int(extraction_range)
	ret_out = open(spacer_mapping_table_url + "_filtered_hits_extracted_faidx_" + "bp_window.fasta", "a")
	next(hit_table)
	for hit in hit_table:
		genome_block = hit[1].split("|")
		genome_id = genome_block[0]
		genome_name = genome_block[1]
		spacer_start = int(hit[8])
		spacer_end = int(hit[9])
		if (spacer_start > spacer_end):
			my_spacer_start = copy.deepcopy(spacer_start)
			spacer_start = spacer_end
			spacer_end = my_spacer_start
		spacer_start = spacer_start - extraction_range
		spacer_end = spacer_end + extraction_range	
		if (spacer_start < 1):
			spacer_start = 1
		# need to label this data with target spacer start and ends!!	
		cmdline = "samtools faidx " + "/g/data/va71/labelled_genomes/" + str(genome_id) + ".fasta " + "\'" + hit[1] + "\'" + ":" + str(spacer_start) + "-" + str(spacer_end)
		cmdline2 = "samtools faidx " + "/g/data/va71/labelled_genomes/" + str(genome_id) + ".fasta " + "\'" + hit[1] + "\'"
		subprocess.run([cmdline],shell=True, stdout=ret_out)
	
	ret_out.close()
	phage_out.close()
	csvfile.close()
	return 0

mapped_genome_retrieval(sys.argv[1])
