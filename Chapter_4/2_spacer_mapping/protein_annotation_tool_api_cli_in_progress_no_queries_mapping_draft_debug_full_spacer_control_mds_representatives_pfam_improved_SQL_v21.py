# protein annotation tool api
import sys
if ('/g/data/va71/crispr_pipeline_annotation' not in sys.path):
	sys.path.append('/g/data/va71/crispr_pipeline_annotation')
# Need to make 5 primary upgrades:
import sqlite3
from Bio import SeqIO
import csv
import os
import subprocess
import pilercr_pos_extractor_annotation
import pilercr_pos_extractor_annotation_full_spacers
import vir_sorter_parser
import hhsuite_parser
import dali_parser
import shutil
import timeit
from time import process_time
from time import time
import fnmatch
import genome_retriever_anno_prot
import rep_cluster_extractor
import genome_identifier_extractor
import table_identifier_extractor
import genome_retriever_target_protein_adder_rep
import protein_hit_of_origin
import phylo_renamer
import pyfasta_but_faster
import hmm_scan_parser
import spacer_dict_generator
from multiprocessing import Pool
import pickle
from annotation_parallelisation_function_SQL8 import annotation_parallelisation
import false_crispr_array_spacers_filtration
import crispr_detect_tabulation
import crisprdetect_inversion
import crispr_mapping_parser_draft_v7_functions_standalone_upgraded
import filtered_arrs_appender_SQL4
import pps_distance_calculator_annotation_corrected_strand_direction
import spacers_2_or_more_2
import genemark_sequence_reformatting
import prodigal_genemark_reconciliation
import genome_representative_clustering
from datetime import date
import pandas
import random
import phage_spacer_table_inversion
import spacer_mapping_parallelisation
import faidx_retriever_genome_block_spacer_mapping_function
import empty_spacer_filterer
import phage_genome_reformatting
import pam_identification
# functions go here!!



# may want to run from command line directly using sys.argv parameters
# taken from: https://stackoverflow.com/questions/1724693/find-a-file-in-python 
def find_pattern(pattern, path):
	result = []
	for root, dirs, files in os.walk(path):
		for name in files:
			if fnmatch.fnmatch(name, pattern):
				result.append(os.path.join(root, name))
		return result		

# relabel all fasta files to include the block_name
def genome_block_relabeller(db_directory_path):
	print("GSTART!!")
	blast_db_seqs = find_pattern('*.fasta', db_directory_path) # this will only work if all files in directory are only fasta!! Need to use find!!
	blast_db_seqs.extend(find_pattern('*.fa', db_directory_path))
	blast_db_seqs.extend(find_pattern('*.fna', db_directory_path))
	# find all fasta files in a directory
	print("loop start!!")
	i=0
	genome_path = blast_db_seqs[0]
	for all_seqs in blast_db_seqs:
		print(all_seqs)
		print(i)
		sequences = SeqIO.parse(all_seqs, "fasta")
		ret_seq = []
		all_seqs = all_seqs.split("/")
		all_seqs = all_seqs[-1]
		my_file = open(db_directory_path + all_seqs + "_labelled.fasta", "w")
		for sequence in sequences:
			this_sequence = sequence
			this_sequence.id =  this_sequence.id # label each db sequence with the block name. This will be important in later steps!!
		#	ret_seq.append(this_sequence)
			SeqIO.write(this_sequence, my_file, "fasta")
		my_file.close()
		print("GEND!!")
		print("good!!")
			# relabel all fasta files to include the block_name
	return genome_path

# helper function to run makeblastdb for db generation
def db_generation_subroutine (db_directory_path):
	# need to format the db fasta file sequences with their respective block names

	genome_path = genome_block_relabeller(db_directory_path) # relabel db files, returning a new file with the labelled prefix attached
	print("db_sub!!")
	subprocess.run(["find " + db_directory_path + " -name *_labelled.fasta | xargs -n 1 -I {} -P 0 makeblastdb -in {} -dbtype nucl"], shell=True) # make BLAST dbs
	print("db_sub_complete!!")
	return genome_path

# run tblastn search to identify matching genomes to query
def db_generation(db_directory_path, file_path):
	print("db_gen_start!!")
	subprocess.run(["find " + db_directory_path + " -name *_labelled.fasta -type f | xargs -n 1 -I {} -P 0 tblastn -query " + file_path + " -db {} -out {}.csv -outfmt 10 -evalue 0.0000001 -max_target_seqs 10000000 -max_hsps 1"], shell=True)
	subprocess.run(["find " + db_directory_path + " -name *_labelled.fasta.csv -type f | xargs -n 1 -I {} -P 1 cat {} >> " + file_path + "_all_hits.csv"], shell=True)
	print("db_gen_end!!")
	return 0


def program_cleaner (): # program to delete all intermediate files generated in order to prevent issues with the program being run twice!!
	return 0

def first_ele (elements):
	return elements[0]

# substitute spaces for hash symbols in protein headers
def hash_adder(input_genome_file_name):
	aa_proteins = list(SeqIO.parse(input_genome_file_name, "fasta"))
	i = 0
	while (i < len(aa_proteins)):
		aa_proteins[i].id = aa_proteins[i].description.replace(" ","#") # getting zero here should be impossible!!
		aa_proteins[i].description = ""
		i += 1
	SeqIO.write(aa_proteins, input_genome_file_name,"fasta")
	return 0

# run protein prediction and reconciliation (if using genemark) on retrieved set of genomes corresponding to each subtype query.
def generate_protein (input_genome_file_name, protein_block_name, rnafold_switch=0,merge_protein=0):
	if (rnafold_switch == 1):
		subprocess.run(["/g/data/va71/crispr_pipeline_annotation/prodigal/prodigal", "-p", "meta", "-i", input_genome_file_name, "-o", input_genome_file_name + "_output_full.txt", "-a", input_genome_file_name + "_aa_raw.fasta", "-d", input_genome_file_name + "_cds_raw.fasta"]) # output should be redirected to genome specific folder
		if (merge_protein == 1):
			subprocess.run(["/g/data/va71/crispr_pipeline_annotation/genemark_reinstall/gms2_linux_64/gms2.pl", "--seq", input_genome_file_name,"--genome-type", "auto", "--out", input_genome_file_name + "_trans_aa_gms2.txt", "--faa", input_genome_file_name + "_trans_aa_gms2.faa", "--fnn", input_genome_file_name + "_trans_aa_gms2.fnn", "--max-iter", "3","--conv-thresh", "0.9"])
			genemark_sequence_reformatting.genemarkS2_reformatting(input_genome_file_name + "_trans_aa_gms2.faa")
			prodigal_genemark_reconciliation.union(input_genome_file_name + "_trans_aa_gms2.faa" + "_geneS2prod_reformatted.fasta", input_genome_file_name + "_aa_raw.fasta", input_genome_file_name + "_aa_raw.fasta") # This will need to be rewritten with additional support to enable fnn sequences to be merged into a single file as well (for RNAfold prediction)
			genemark_sequence_reformatting.genemarkS2_reformatting(input_genome_file_name + "_trans_aa_gms2.fnn")
			prodigal_genemark_reconciliation.union(input_genome_file_name + "_trans_aa_gms2.fnn" + "_geneS2prod_reformatted.fasta", input_genome_file_name + "_cds_raw.fasta", input_genome_file_name + "_CDS_raw.fnn")
		subprocess.run(["/g/data/va71/crispr_pipeline_annotation/seqkit", "faidx", input_genome_file_name + "_CDS_raw.fnn"])
	else:		
		subprocess.run(["/g/data/va71/crispr_pipeline_annotation/prodigal/prodigal", "-p", "meta", "-i", input_genome_file_name, "-o", input_genome_file_name + "_output_full.txt", "-a", input_genome_file_name + "_aa_raw.fasta" ]) # output should be redirected to genome specific folder
		if (merge_protein == 1):
				subprocess.run(["/g/data/va71/crispr_pipeline_annotation/annotation_database/annotation_upgrades/protein_reconciliation/genemark_reinstall/gms2_linux_64/gms2.pl", "--seq", input_genome_file_name,"--genome-type", "auto", "--out", input_genome_file_name + "_trans_aa_gms2.txt", "--faa", input_genome_file_name + "_trans_aa_gms2.faa", "--max-iter", "3","--conv-thresh", "0.9"])
				genemark_sequence_reformatting.genemarkS2_reformatting(input_genome_file_name + "_trans_aa_gms2.faa")
				prodigal_genemark_reconciliation.union(input_genome_file_name + "_trans_aa_gms2.faa" + "_geneS2prod_reformatted.fasta", input_genome_file_name + "_aa_raw.fasta", input_genome_file_name + "_aa_raw.fasta") # This will hopefully override the original protein file
	# This is where samtools is called. Need to add a function on input_protein_file to replace the whitespace with an allowed character.
	subprocess.run(["/g/data/va71/crispr_pipeline_annotation/seqkit", "faidx", input_genome_file_name + "_aa_raw.fasta"])	
	subprocess.run(["makeblastdb", "-in", input_genome_file_name + "_aa_raw.fasta", "-dbtype", "prot"])
	subprocess.run(["cat " + input_genome_file_name + "_aa_raw.fasta >> " + protein_block_name],shell=True)	# This may result in some protein duplication, however samtools does not duplicate the indexes.
	subprocess.run(["/g/data/va71/crispr_pipeline_annotation/seqkit", "rmdup",protein_block_name])
	subprocess.run(["/g/data/va71/crispr_pipeline_annotation/seqkit", "faidx","--update-faidx", protein_block_name])
	return 0

# script to reverse the order of the entries post reconcilation (if reverse is selected)
def table_reverse(input_table_url):
	with open(input_table_url) as csvfile:
		hit_table = list(csv.reader(csvfile))
	table_dict = {}
	for hit in hit_table:
		if hit[0] not in table_dict:
			table_dict[hit[0]] = [hit]
		else:
			table_dict[hit[0]].append(hit)
	hit_table = table_dict.values()
	i = 0
	while (i < len(hit_table)):
		if(hit_table[i][0][1] == "Reverse"):
			hit_table[i].reverse()
		i += 1
	ret_out = open(input_table_url, "w")
	spam_writer = csv.writer(ret_out)
	for genome in hit_table:
		for row in genome:
			spam_writer.writerow(row)
	ret_out.close()
	return 0		

# create genomes table from a set of genomes without an input query.
def filler_table_generation (my_genomes, genome_url, params, sql_db_connect):
	# need to lookup GOLD db to retrieve the matching Analysis ID if an NCBI id is supplied
	con = sql_db_connect
	cur = con.cursor()
	ret_list = []
	header = ["RUN", "Seed_id", "Genome_id", "perc_identity", "length", "mismatch", "gapopen", "query_start", "query_end", "match_start", "match_end","evalue","bitscore", "my_date", "params","Genome_id_prefix"]
	ret_list.append(header)
	out_ret = open(genome_url + "_all_hits.csv", "w")
	spam_writer = csv.writer(out_ret)
	spam_writer.writerow(header)

	for genome in my_genomes:
		# This only works if all genomes are JGI derived!! Need to create function to lookup genome_id pseudonyms 
		genome_id_prefix = genome.id.split("_")[0]
		if (genome_id_prefix[:1] == "Ga"):
			pass 
		else:
			genome_id_prefix = None

		row = ["Not assigned",None,genome.id, None, None,None,None, None,None,None, None,None,None,date.today(),params, genome_id_prefix]
		spam_writer.writerow(row)
		ret_list.append(row)	
	con.close()	
	return 0

# generate table compatible with SQL. Note: this function was not used
def genome_hit_sql_table_generation(hit_table, genome_url, params, sql_db_connect):
	ret_list = []
	header = ["RUN","Seed_id", "Genome_id", "perc_identity", "length", "mismatch", "gapopen", "query_start", "query_end", "match_start", "match_end","evalue","bitscore", "my_date", "params","Genome_id_prefix","Genome_path"]
	ret_list.append(header)
	out_ret = open(genome_url, "w")
	spam_writer = csv.writer(out_ret)
	spam_writer.writerow(header)
	print(header)
	for hit in hit_table:
		genome_id_prefix = hit[1].split("_") [0]
		if (genome_id_prefix[:1] == "Ga"):
			row = ["Not_assigned"] + hit + [date.today(), params] + [genome_id_prefix]
		else:
			row = ["Not_assigned"] + hit + [date.today(), params] + [None] # These may need to be null values	
		spam_writer.writerow(row)
		ret_list.append(row)	
	return 0	

# Retrieve genome sequences by identifier.
def genome_hit_lookup (hit_table, db_directory_path, genome_file_url, block, genome_row_no=2): # b?
	row_dict = {}
	if (os.path.isfile(genome_file_url)):
		os.remove(genome_file_url)
	genome_file = open(genome_file_url, "a")
	for row in hit_table:

		genome_row = row[genome_row_no]
		genome_row = genome_row.split("|") 
		true_gene_id = genome_row[0]

		genome_row = block
		if (genome_row not in row_dict):
			my_id = {true_gene_id}
			row_dict[genome_row] = my_id
		else:
			row_dict[genome_row].add(true_gene_id) 
	for row in row_dict.items():
		genome_row = row[0]
		print(genome_row)
		block = SeqIO.parse(genome_row, "fasta")
		for gene in block:
			gene_key = gene.id.split("|")
			if (len(gene_key) == 1):
				gene_key = gene_key[0]
			else:
				gene_key = gene_key[1]	
			if (gene_key in row[1]):
				my_genome = gene
				my_genome.id = gene_key
				my_genome.description = ""  # my_genome.description # + "|" + row[0] # now sequence should be annotated with block + genome_id + sequence_id
				SeqIO.write(my_genome, genome_file, "fasta") # This needs to be written to file then indexed with samtools
	genome_file.close()
	return 0

#START!!
# This is the main code for performing spacer mapping
# input_flags (full set of options was not implemented):
# -i input (either csv_file or input sequence in fasta format)
# -t input_type (query/table)
# -p use phmmer
# -b use hhblits
# -s use hhsearch
# -a use alphafold + DALI
# should supply a number of potential tokens as well as having default values. Tokens are represented with an '-' specifier.

# INPUT: 1. representative ORF/ or other CRISPR-associated protein sequence (in FASTA format)
# 		 i.e. cas13a.fasta (located in "queries/" folder)
# 		 2. file path to folder containing 10TB genome blocks, along with  files to run BLAST, for each block
#        i.e. /g/data/va71/labelled_genomes/
#		 3. file path to folder called "genomes/" containing the file with DNA extracted 20kb upstream and downstream of CRISPR-arrays
#		 i.e. "cas13a/genomes/" containing "spliced_debug_corrected_cleaned_rerun_20kb_local_windows_combined.fasta" and associated files for running BLAST
#		 4. Additionally, in order to run correctly, the output from running "" must be present in the "queries/" and "genomes/" directories
# OUTPUT: 1. Table containing hits to mapped spacers. This table has been filtered to exclude hits to self-CRISPR arrays, and hits to other arrays encodding the same spacers
#		  2. FASTA file containing a list of "host" encoding contigs
#		  3. FASTA file containing spacer mapped target contigs (in FASTA format)
#		  4. Preliminary table of PPS-spacer distances (not used)
#		  5. Preliminary PAM concensous predictions for each subtype (not used)


# SHELL: see annotation_cas13a_working_spacer_mapping1.sh. located in "2_spacer_mapping/running_scripts/"

i = 1

db_directory_path = "genomes/"
db_directory_switch = 0
phmmer_switch = 0
hhsearch_switch = 0
hhblits_switch = 0
dali_switch = 0
spacer_mapping_switch = 0
block_directory = " /g/data/va71/labelled_genomes/"
file_path = ""
bypass_switch = 0
target_pos_switch = 0
protein_query_switch = 0
criscasfinder_switch = 0
pfam_switch = 0
large_dataset = 0
phylogeny_switch = 0
spacer_only_switch = 0
min_seq_id = 0.4
tony_mapping_switch = 0
mapping_skip_switch = 0
crispr_orientation = 0
crispr_detect = 0
merge_protein = 0
phage_genomes = 0
pam_anno = 0
gold_switch = 0
phage_protein_query_switch = 0
protein_generation_switch = 0
cores = 12
dont_initialise = 0
spacer_generation_bypass_switch = 0
eggnog_switch=0
vmatch_switch=0
virsorter_switch=0
rnafold_switch=0
padlocplus_switch_novel=0



# may want to add switches for phmmer/hhsearch/hhblits database types below
print("welcome!!")

# may also want to parse output types and other modifying parameters. For now these will be omitted.
while (i < len (sys.argv)):
	match sys.argv[i]:
		case '-i':
			if (1 + i > len(sys.argv) or sys.argv[i + 1] in ("-i", "-t", "-d", "-m", "-n") ):
				print("no input argument given")
				exit()
			else:
				file_path = sys.argv[i + 1]	# protein sequence in FASTA format
			i += 2	
		case '-t':
			if (1 + i > len(sys.argv) or sys.argv[i + 1] in ("-i", "-t", "-d", "-m", "-n")):
				print("no_input type argument given")
				exit()

			else:
				input_type = sys.argv[1 + i]
				if (input_type == "query"): # change to quer(ies)?
					db_directory_switch = 1

				elif (input_type == "table"):
					db_directory_switch = 2
				else:
					print("incorrect input_specified!!")
					exit()
			i += 2	
		case '-m':
			if (1 + i > len(sys.argv) or sys.argv[i + 1] in ("-i", "-t", "-d", "-m", "-n")):
				print("no block directory given!!")
				exit()
			else:
				spacer_mapping_switch = 1
				block_directory = sys.argv[i + 1] 	
			i += 2		
		case '-d':
			if (1 + i > len(sys.argv) or sys.argv[i + 1] in ("-i", "-t", "-d", "-m", "-n")):
				print("no db directory path given!!")
			else:
				db_directory_path = sys.argv[i + 1]	
			i += 2
		case '-n':
			if (1 + i) > len(sys.argv[i + 1] in ("-i", "-t", "-d", "-m", "-n")):
				print("no min_seq_id_value_given!!")
				exit()
			else:
				min_seq_id = int(sys.argv[i + 1])
			i += 2					
		case '-x':
			bypass_switch = 1	
			i += 1
		case '-p':
			phmmer_switch = 1
			i += 1
		case '-s':
			hhsearch_switch	= 1
			i += 1
		case '-b':
			hhblits_switch = 1
			i += 1	
		case '-l':
			dali_switch = 1
			i += 1	
		case '-r':
			target_pos_switch = 1
			i += 1
		case '-q':
			protein_query_switch = 1
			i += 1
		case '-f': # annotate using pfam database
			pfam_switch = 1
			i += 1
		case '-c': # annotate using CRISPRCasFinder profiles
			criscasfinder_switch = 1
			i += 1 	
		case '-L':
			large_dataset = 1
			i += 1	
		case '-p':
			phylogeny_switch = 1
			i += 1	
		case '-e':
			spacer_only_switch = 1	
			i += 1
		case '-tony':
			tony_mapping_switch = 1
			i += 1	
		case '-skip':
			mapping_skip_switch = 1	
			i += 1
		case '-orientation':
			crispr_orientation = 1	
			i += 1
		case '-crispr_detect':
			crispr_detect = 1
			i += 1	
		case '-merge_prot':
			merge_protein = 1	
			i += 1
		case '-phage_genomes':
			phage_genomes = 1
			i += 1
		case '-PAM_anno':
			pam_anno = 1
			i += 1	
		case '-ignore_initialisation':
			dont_initialise = 1
			i += 1
				
		case '-phage_rep':
			phage_protein_query_switch = 1	
			i += 1
		case '-prot_case_generation':
			protein_generation_switch = 1
			i += 1	
		case '-add_gold':
			gold_switch =  1
			i += 1	
		case '-cores':
			if ((1 + i) > len(sys.argv[i + 1]) in ("-i", "-t", "-d", "-m", "-n")):
				print("no min_seq_id_value_given!!")
				exit()
			else:
				cores = int(sys.argv[1+i])
			i += 2
		case '-spacer_gen_bypass':
			spacer_generation_bypass_switch = 1
			i += 1
		case '-eggnog':
			eggnog_switch = 1
			i += 1
		case '-vmatch_switch':
			vmatch_switch = 1
			i += 1
		case '-virsorter_switch':
			virsorter_switch = 1
			i += 1
		case '-rnafold_switch':
			rnafold_switch = 1
			i += 1
		case '-padlocplus_switch_novel':
			padlocplus_switch_novel = 1
			i += 1		

eggnog_switch=0
vmatch_switch=0
virsorter_switch=0
rnafold_switch=0
padlocplus_switch_novel=0			


print("checkpoint 0")

sql_db_connect = sqlite3.connect("/g/data/va71/crispr_pipeline_annotation/" + "crispr-phage_interactome.sql")
genome_hits_block = "/g/data/va71/crispr_pipeline_annotation/" + "genome_annotation_block.fasta" 
phage_hits_block = "/g/data/va71/crispr_pipeline_annotation/" + "phage_annotation_block.fasta"
protein_hits_block = "/g/data/va71/crispr_pipeline_annotation/" + "protein_annotation_block.fasta"
phage_protein_hits_block = "/g/data/va71/crispr_pipeline_annotation/" + "phage_protein_annotation_block.fasta"
whole_phage_block_url = "/g/data/va71/crispr_pipeline_annotation/" + "whole_phage_genome_block.fasta"
cur = sql_db_connect.cursor()

golddata_url = "/g/data/va71/crispr_pipeline_annotation/golddata/golddata_7_6_2023.xlsx"
if (spacer_generation_bypass_switch == 0):
	if (gold_switch == 1):
		gold_study = pandas.read_excel(golddata_url,sheet_name='Study')
		gold_biosample = pandas.read_excel(golddata_url,sheet_name='Biosample')
		gold_Organism = pandas.read_excel(golddata_url,sheet_name='Organism')
		gold_Sequencing_Project = pandas.read_excel(golddata_url,sheet_name='Sequencing Project')
		gold_Analysis_project = pandas.read_excel(golddata_url,sheet_name="Analysis Project")


	# Intialise SQL db (functionality was not used in final runs)
		try:
			gold_study.to_sql("GOLD_STUDY",sql_db_connect,if_exists='fail',index=False)
			cur.execute("CREATE INDEX GOLD_STUDY_ID_INDEX ON GOLD_STUDY (`STUDY GOLD ID`)")
		except ValueError:	
			pass 
		try:
			gold_biosample.to_sql("GOLD_BIOSAMPLE",sql_db_connect,if_exists='fail',index=False)
			cur.execute("CREATE INDEX BIOSAMPLE_GOLD_ID_INDEX ON GOLD_BIOSAMPLE (`BIOSAMPLE GOLD ID`)")
			cur.execute("CREATE INDEX GOLD_BIOSAMPLE_NCBI_TAX_ID ON GOLD_BIOSAMPLE (`BIOSAMPLE NCBI TAX ID`)")
		except ValueError:
			pass
		try:
			gold_Organism.to_sql("GOLD_ORGANISM", sql_db_connect,if_exists='fail',index=False)
			cur.execute("CREATE INDEX ORGANISM_GOLD_ID_INDEX ON GOLD_ORGANISM (`ORGANISM GOLD ID`)")
			cur.execute("CREATE INDEX ORGANISM_NCBI_TAX_ID_INDEX ON GOLD_ORGANISM (`ORGANISM NCBI TAX ID`)")
		except ValueError:
			pass 
		try:
			gold_Sequencing_Project.to_sql("GOLD_SEQUENCING_PROJECT", sql_db_connect,if_exists='fail',index=False)
			cur.execute("CREATE INDEX PROJECT_GOLD_ID_INDEX ON GOLD_SEQUENCING_PROJECT (`PROJECT GOLD ID`)")
			cur.execute("CREATE INDEX PMO_PROJECT_ID_INDEX ON GOLD_SEQUENCING_PROJECT (`PMO PROJECT ID`)")
			cur.execute("CREATE INDEX NCBI_BIOPROJECT_ACCESSION_INDEX ON GOLD_SEQUENCING_PROJECT (`NCBI BIOPROJECT ACCESSION`)")
			cur.execute("CREATE INDEX NCBI_BIOSAMPLE_ACCESSION_INDEX ON GOLD_SEQUENCING_PROJECT (`NCBI BIOSAMPLE ACCESSION`)")
			cur.execute("CREATE INDEX SRA_EXPERIMENT_IDS_INDEX ON GOLD_SEQUENCING_PROJECT (`SRA EXPERIMENT IDS`)")
			cur.execute("CREATE INDEX PROJECT_LEGACY_GOLD_ID_INDEX ON GOLD_SEQUENCING_PROJECT (`PROJECT LEGACY GOLD ID`)")
			cur.execute("CREATE INDEX STUDY_GOLD_ID_INDEX ON GOLD_SEQUENCING_PROJECT (`STUDY GOLD ID`)")

		except ValueError:
			pass 
		try: 
			gold_Analysis_project.to_sql("GOLD_ANALYSIS_PROJECT",sql_db_connect,if_exists='fail',index=False)
			cur.execute("CREATE INDEX AP_GOLD_ID_INDEX ON GOLD_ANALYSIS_PROJECT (`AP GOLD ID`)")
			cur.execute("CREATE INDEX AP_ITS_AP_ID_INDEX ON GOLD_ANALYSIS_PROJECT (`AP ITS AP ID`)")
			cur.execute("CREATE INDEX AP_IMG_TAXON_ID_INDEX ON GOLD_ANALYSIS_PROJECT (`AP IMG TAXON ID`)")
			cur.execute("CREATE INDEX AP_GENBANK_INDEX ON GOLD_ANALYSIS_PROJECT (`AP GENBANK`)")
			cur.execute("CREATE INDEX AP_NCBI_SRA_IDS_INDEX ON GOLD_ANALYSIS_PROJECT (`AP NCBI SRA IDS`)")
			cur.execute("CREATE INDEX AP_STUDY_GOLD_ID_INDEX ON GOLD_ANALYSIS_PROJECT (`AP STUDY GOLD ID`)")
			cur.execute("CREATE INDEX AP_ORGANISM_GOLD_ID ON GOLD_ANALYSIS_PROJECT (`AP ORGANISM GOLD ID`)")
			cur.execute("CREATE INDEX AP_PROJECT_GOLD_IDS ON GOLD_ANALYSIS_PROJECT (`AP PROJECT GOLD IDS`)")
		except ValueError:
			pass			


		sql_db_connect.commit()

	print("Opened_database_successfully")
	# if no input file is given, the program will default to using all the genomes as input.
	if '-i' not in sys.argv:
		db_directory_switch = 0
	print(db_directory_switch)
	# generate db to retrieve genomes via BLAST search
	if (dont_initialise == 0):
		if (os.path.isdir(db_directory_path)):
			if (db_directory_switch == 1):
				if not (os.path.isfile(db_directory_path + "_labelled.fasta")): # This is NOT the filename. Hard to retrieve without a grep-like command
					genome_path = db_generation_subroutine(db_directory_path) # relabel files and compile blast db using genomes as input
					db_generation(db_directory_path, file_path)	# run blast to generate result file in csv format
				else:
					print("error in db generation")
					exit()
			elif (db_directory_switch == 2):
				if not (os.path.isfile(db_directory_path + "_labelled.fasta")):
					genome_path = db_generation_subroutine(db_directory_path)

			elif(db_directory_switch == 0):
				print("good!!")
				genome_path =genome_block_relabeller(db_directory_path)

			else:
				print("error: incorrect subroutine_option!!")		
		else:
			print("error:wrong BLAST db directory path supplied!!")	
			print(db_directory_path)
			exit()


	print("checkpoint 1!!")

	# retrieve the corresponding genomes using query-based BLAST search or a table of genome sequence IDs.
	print("db_directory_switch:")
	print(db_directory_switch)
	if (db_directory_switch != 0):
		if (input_type == "query"): # here we want to take the table generated from db_generation
			mod_file_path = file_path # .split('/')
			a = mod_file_path + "_all_hits.csv" #  This should point to the generated table. May need a '/' between db_directory_path and file path for this to work!!
		elif (input_type == "table"): # here we want to take the input table, which has been labelled!!
			a = mod_file_path # this will be a specified table.
		else: 
			print("problem with the input character!!")



		with open (a, "r") as csvfile:
			hit_table = list(csv.reader(csvfile))

		b = a + "_genomes.fasta" 
		genome_hit_lookup(hit_table, db_directory_path, b, genome_path, 1)
		genome_file = open(b, "r")

		# extract the block name from each table row
		previous_genome_row = ""
		# if samtools file exists then delete file. Else samtools and index the genomes
		if (os.path.isfile(genome_hits_block + ".fai")):
			# need to only add genomes which are not already in the data block. Extract the ids from the samtools file
			with open (genome_hits_block + ".fai", "r") as csvfile:
				fai_table = list(csv.reader(csvfile,delimiter='\t'))
			subprocess.run(["rm " + genome_hits_block + ".fai"],shell=True)	
			fai_headers = set()
			for row in fai_table:
				if (row[0] not in fai_headers):
					fai_headers.add(row[0])
				else: 
					print("Error!! Should not have duplicate fai headers in one job!!")
			new_genomes = SeqIO.parse(genome_file, "fasta")
			existing_block = open(genome_hits_block, "a")
			for genome in new_genomes:
				print(genome.id)
				if (genome.id not in fai_headers):
					SeqIO.write(genome, existing_block, "fasta")
			existing_block.close()


		else: # start with a set of input genomes, as opposed to a query sequence	
			subprocess.run(["cat " + b + " >> " + genome_hits_block],shell=True)
		subprocess.run(["/g/data/va71/crispr_pipeline_annotation/seqkit rmdup " + genome_hits_block],shell=True)
		subprocess.run(["/g/data/va71/crispr_pipeline_annotation/seqkit faidx " + "--update-faidx " + genome_hits_block],shell=True)
		print("Done!!")	
		# if option -q is used then prodigal should be used to predict proteins arising from tblastn hits.
		# blastp should then identify these proteins.
		# these proteins should then be retrieved!!
		# these should also be a representative set from mmseqs clustering, so that redundant proteins are removed.
		
		if (protein_query_switch == 1):
			representative_table = genome_representative_clustering.rep_cluster(b, min_seq_id, db_directory_path, a, merge_protein, large_dataset, phylogeny_switch, sys.argv)
			a = representative_table
			with open (a, "r") as csvfile:
				hit_table = list(csv.reader(csvfile))
		
		# generate SQL table (vestigal function)
		genome_hit_sql_table_generation(hit_table,a, "".join(sys.argv), sql_db_connect) # need to reach here to bypass creating new genome blocks. The reason is because this function modifies the hit_table!!!
		
		# RUN this code as a standalone. Attempt to integrate the prefix ids. Should not need to rest!!!!
		export_table = pandas.read_csv(a,header=0)
		genome_prefix_ids = {}
		genome_type = {}
		genome_is = "HOST"
		for entry in export_table['Genome_id']:
			my_entry = entry
			my_entry = my_entry.split("_") 	
			if (len(my_entry) > 1):
				my_entry = my_entry[0]
			#	print(my_entry)
				if 'Genome_id_prefix' not in genome_prefix_ids:
					genome_prefix_ids['Genome_id_prefix'] = [my_entry]
					genome_type['Genome_type'] = [genome_is]
				else:
					genome_prefix_ids['Genome_id_prefix'].append(my_entry)
					genome_type['Genome_type'].append(genome_is)
			else:
				if 'Genome_id_prefix' not in genome_prefix_ids:
					genome_prefix_ids['Genome_id_prefix'] = [None]
					genome_type['Genome_type'] = [None]
				else:
					genome_prefix_ids['Genome_id_prefix'].append(None)
					genome_type['Genome_type'].append(None)
							
		export_table = export_table.assign(Genome_path = genome_path)	
		export_table = export_table.assign(Genome_id_prefix = list(genome_prefix_ids.values()) [0]) # These don't seem to work!
		export_table = export_table.assign(Genome_type = list(genome_type.values()) [0]) # These don't seem to work!
		
		cur.execute("PRAGMA foreign_keys = ON;")
		cur.execute('CREATE TABLE IF NOT EXISTS GENOMES (RUN STRING, Seed_id STRING, Genome_id STRING PRIMARY KEY, perc_identity STRING, length STRING, mismatch STRING, gapopen STRING, query_start STRING, query_end STRING, match_start STRING, match_end STRING, evalue STRING, bitscore STRING, my_date STRING, params STRING, Genome_type STRING, Genome_id_prefix TEXT, Genome_path TEXT)') # FOREIGN KEY (Genome_id_prefix) REFERENCES GOLD_ANALYSIS_PROJECT (`AP GOLD ID`) see if this works without the foreign key. Might not be feasible to make mapping 1 to 1 due to NCBI IDs.
		sql_db_connect.commit()
		# generate a unique key or keep trying
		while (True):
			index = random.randint(1000000000000000000, 10000000000000000000)
			index = str(index)
		#	index = random.choices(L, k = 1) [0]
			res = cur.execute('SELECT * FROM GENOMES WHERE RUN LIKE ?', (index,))
			row = res.fetchone()
			if (row != None): # 
				print("random_generation!!") 
			else:
				break
					

	else: # need an if statement to delete this directory prior to running the program. Ditto for BLAST compiled db!!

		subprocess.run(["find " + db_directory_path + " -name *_labelled.fasta -type f | xargs -n 1 -I {} -P 1 cat {} >> " + db_directory_path  + "genome_only_run_master_file.fasta"], shell=True) # first concatenate the genomes
		a = db_directory_path + "genome_only_run_master_file.fasta"
		b = a + "_genomes.fasta"
		my_genomes = SeqIO.parse(a, "fasta")
		filler_table_generation(my_genomes, a, "".join(sys.argv), sql_db_connect)
		export_table = pandas.read_csv(a + "_all_hits.csv")
		genome_prefix_ids = {}
		genome_type = {}
		genome_is = "HOST"
		for entry in export_table['Genome_id']:
			my_entry = entry.split("_") [0]
			if 'Genome_id_prefix' not in genome_prefix_ids:
				genome_type['Genome_type'] = [genome_is]
			else:
				genome_type.append(genome_type)
			
		export_table['Genome_type'] = list(genome_type.values()) [0]

		cur.execute('CREATE TABLE IF NOT EXISTS GENOMES (RUN STRING, Seed_id STRING, Genome_id STRING PRIMARY KEY, perc_identity STRING, length STRING, mismatch STRING, gapopen STRING, query_start STRING, query_end STRING, match_start STRING, match_end STRING, evalue STRING, bitscore STRING, my_date STRING, params STRING, Genome_type STRING, Genome_id_prefix STRING, Genome_path TEXT, FOREIGN KEY (Genome_id_prefix) REFERENCES GOLD_ANALYSIS_PROJECT (`AP GOLD ID`))')
		cur.execute('CREATE INDEX IF NOT EXISTS GENOME_PREFIX_INDEX ON GENOMES(Genome_id_prefix)')
		cur.commit()
		export_table.to_sql('GENOMES',sql_db_connect, if_exists='append', index=True)


		SeqIO.write(my_genomes, b, "fasta")

	sql_db_connect.close()
	print("checkpoint 2!!")

	# write prodigal + genemark orf translation of target genomes here!!
	# also need to index with samtools!!
	if (protein_generation_switch == 1):
		input_genome_file_name = b
		generate_protein(input_genome_file_name,protein_hits_block, rnafold_switch, merge_protein)
	print("checkpoint 3!!")

	if (spacer_mapping_switch == 1):
		print("checkpoint 4!!")
		if (os.path.isdir(db_directory_path + "spacer_distribution_analysis/")):
				pass 
		else:	
			os.mkdir(db_directory_path + "spacer_distribution_analysis/")
		b_basename = b.split("/")
		b_basename = b_basename[-1]

		# need to seperate split spacers into smaller files
		output_dir = db_directory_path + "spacer_distribution_analysis/" + b_basename
		partition_size = os.path.getsize(b) // cores
		partitions = pyfasta_but_faster.rename(b,output_dir,  1966080) # was 1966080
		m = 0
		# may be able to parallelise this loop to save KSU!!
		partition_cores = 8
		protopool = Pool(int(cores // partition_cores))
		
	# This section will probably need to be parallelised to compensate for the dramatically increased run time due to soo many prediction tools running.
		while (m <= partitions):
			partition_dir = output_dir + "_partition_" + str(m) + ".fasta"
			protopool.apply_async(spacer_mapping_parallelisation.parallel_spacer_partition, (output_dir,partition_dir, db_directory_path,b_basename,crispr_detect, crispr_orientation, partition_cores))
			m += 1
		protopool.close()
		protopool.join()
else: # execute this if CRISPR-arrays have already been predicted and extracted.

	b_basename = file_path + "_all_hits.csv" + "_genomes.fasta"
	b_basename = b_basename.split("/") [-1]
	# Need to put in a bypass to jump to this line of code!!!!!!

# predict CRISPR-arrays and merge predictions (run seperately of the main workflow intially)
if (crispr_detect == 1 and crispr_orientation == 1):
	pilercr_pos_extractor_annotation_full_spacers.merged_detectcrtpilercr_table_2_fasta(db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs" + ".lst" + "_reconciled_full_arr_positions.csv")
else:	
	pilercr_pos_extractor_annotation_full_spacers.spacer_table_2_fasta(db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs" + ".lst" + "_full_real_arr_positions.csv")
	 
if (tony_mapping_switch == 1):
	offset = 10
	if (crispr_detect == 1 and crispr_orientation == 1):
		pilercr_pos_extractor_annotation_full_spacers.merged_detectcrtpilercr_table_2_fasta_with_dr_offset_all_cases_but_better(db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs" + ".lst" + "_reconciled_full_arr_positions.csv", offset)
	else:	
		pilercr_pos_extractor_annotation_full_spacers.spacer_table_2_fasta_with_dr_offset_all_cases(db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs" + ".lst" + "_full_real_arr_positions_tony.csv", offset)
	
#	append the header
with open(db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs" + ".lst" + "_reconciled_full_arr_positions.csv","r") as csvfile:
	hit_table = list(csv.reader(csvfile))
out_file = open(db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs" + ".lst" + "_reconciled_full_arr_positions.csv","w")
spam_writer = csv.writer(out_file)
spam_writer.writerow(["Genome_id","orientation","orientation_score", "orientation_confidence", "questionable_array", "array_score","CRISPR-start","CRISPR-end", "repeat_start", "repeat_end","spacer_start","spacer_end","dr_repeat_original", "dr_repeat_concensous", "spacer", "Array_tool" ])
for hit in hit_table:
	spam_writer.writerow(hit)
out_file.close()	
# Parameters for spacer mapping!!
# These parameters should be related to the dr offset score!!!
perc_identity = 0.90
dr_perc_identity = 0.90
query_cover = 95
dr_query_hsp_cover = 95
formatting = "10 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen"
if (mapping_skip_switch == 0):
	if (crispr_detect == 1 and crispr_orientation == 1):
		input_spacer_fasta = db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs" + ".lst" + "_reconciled_full_arr_positions.csv" + "_reconciled_spacers.fasta"
	else:
		input_spacer_fasta = db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs" + ".lst" + "_full_real_arr_positions.csv" + "_spacers.fasta"
	# filter spacers which are misidentified and actually concatenated elements of entire arrays
	false_crispr_array_spacers_filtration.filtration(input_spacer_fasta)
	input_spacer_fasta = input_spacer_fasta + "_arrs_filtered.fasta"
	if (tony_mapping_switch == 1):
		if (crispr_detect == 1 and crispr_orientation == 1):
			input_dr_spacer_fasta = db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs" + ".lst" + "_reconciled_full_arr_positions.csv" + "_drs_" + str(offset) + "_reconciled_spacers.fasta"
		else:
			input_dr_spacer_fasta = db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs" + ".lst" + "_full_real_arr_positions_tony.csv" + "_drs_" + str(offset) + "_spacers.fasta"
		# filter spacers which are misidentified and actually concatenated elements of entire arrays
		false_crispr_array_spacers_filtration.filtration(input_dr_spacer_fasta)
		input_dr_spacer_fasta = input_dr_spacer_fasta + "_arrs_filtered.fasta"
		print("Good!!")
		subprocess.run(["/g/data/va71/crispr_pipeline_annotation/pipe_line_source_files/nested_mpirun_script_no_qsub_flock_optimised_test_script.sh", input_dr_spacer_fasta, block_directory, str(dr_perc_identity), str(cores), "/g/data/va71/crispr_pipeline_annotation/pipe_line_source_files/", db_directory_path, "co_occurrance_bug_diagnostics_results2/", str(dr_query_hsp_cover), formatting]) # what input args are required??
	print("Pipeline parameters:")
	print(input_spacer_fasta)
	print (block_directory)
	print (perc_identity)
	print(cores)
	# subproccess call to script which maps spacers barcoded with direct repeat handles
	subprocess.run(["/g/data/va71/crispr_pipeline_annotation/pipe_line_source_files/nested_mpirun_script_no_qsub_flock_optimised_test_script.sh", input_spacer_fasta, block_directory, str(perc_identity), str(cores), "/g/data/va71/crispr_pipeline_annotation/pipe_line_source_files/", db_directory_path, "co_occurrance_bug_diagnostics_results/", str(query_cover), formatting ]) # what input args are required??
		
if (tony_mapping_switch == 1):
	# exclude mapped hits to crispr arrays (or degraded array homologs)
	subprocess.run(["/g/data/va71/crispr_pipeline_annotation/pipe_line_source_files/block_tony_spacer_filtration5.sh", db_directory_path])	
	print ("Up to penultimate step!!:")
	print("folder dir:")
	print(db_directory_path)
	print("genomes_dir:")
	print(b_basename)

	subprocess.run(["/g/data/va71/crispr_pipeline_annotation/result_cons_script.sh", db_directory_path, b_basename])
# Write merging with CRISPR_Info here!!
crispr_array_table = db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs" + ".lst" + "_reconciled_full_arr_positions.csv"
array_table = pandas.read_csv(crispr_array_table)
i = 0
run_column = {"RUN":[]}
while ( i < len(array_table["Genome_id"])):	
	run_column["RUN"].append(index) # This variable is randomly generated a long way upstream in the code. While an index can be generated to simulate this. It won't be the same as running from the beginning
	i += 1

array_table["RUN"] = run_column["RUN"]
# In this table the primary keys will be the spacers/ spacer-genome-id pairs as these are unique
# This table needs to be interfaced with using foreign keys. But the really important table for investigations is definitely the Filtered spacer mapping table!!
sql_db_connect = sqlite3.connect("/g/data/va71/crispr_pipeline_annotation/" + "crispr-phage_interactome.sql")
cur = sql_db_connect.cursor()
cur.execute("CREATE TABLE IF NOT EXISTS RAW_ARRAYS (RUN STRING, Genome_id STRING, orientation STRING, orientation_score STRING, orientation_confidence STRING, questionable_array STRING, array_score STRING, `CRISPR-start` STRING,`CRISPR-end` STRING,repeat_start STRING, repeat_end STRING, spacer_start STRING, spacer_end STRING, dr_repeat_original STRING, dr_repeat_concensous STRING, spacer STRING, Array_tool STRING, FOREIGN KEY (Genome_id) REFERENCES GENOMES (Genome_id) )")
cur.execute("CREATE INDEX IF NOT EXISTS GENOME_ID_INDEX ON RAW_ARRAYS (`Genome_id`)")
sql_db_connect.commit()
array_table.to_sql('RAW_ARRAYS',sql_db_connect, if_exists='append', index=False) # May need to make a table tabulating all the run numbers, seeds and dates run.
array_table.to_csv(crispr_array_table, index=False)
empty_spacer_filterer.remove_blanks(db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs" + ".lst" + "_reconciled_full_arr_positions.csv")

# expand and number spacers in each concensus CRISPR array
spacer_hitmap_master_table = filtered_arrs_appender_SQL4.spacer_hits_appender(db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered.csv", db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs" + ".lst" + "_reconciled_full_arr_positions.csv_no_blanks.csv", db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_hitmap.csv")
hitmap_table_url = db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_hitmap.csv"
hit_table = pandas.read_csv(hitmap_table_url)
	
# Need to identify whether hits to phage genomes have already occurred for seperate runs. Are these genomes 

# CODE to transfer table to SQL database. Note, this functionality was not used other than to check subtype sequences did not overlap.
cur.execute("CREATE TABLE IF NOT EXISTS SPACER_HITMAP (Spacer_id TEXT, Phage_id TEXT, Perc_id FLOAT, Length INT, Mismatches TEXT, Gapopen TEXT, query_start TEXT, query_end TEXT, Mapped_start_site TEXT, Mapped_end_site TEXT, evalue TEXT, bitscore TEXT, Genome_id TEXT, orientation TEXT, orientation_score TEXT, orientation_confidence TEXT, questionable_array TEXT, array_score TEXT, `CRISPR-start` TEXT,`CRISPR-end` TEXT,repeat_start TEXT, repeat_end TEXT, spacer_start TEXT, spacer_end TEXT, dr_repeat_original TEXT, dr_repeat_concensous TEXT, spacer TEXT, Array_tool TEXT, RUN TEXT, array_number TEXT, spacer_number TEXT, mapped_genomes TEXT, FOREIGN KEY (Genome_id) REFERENCES GENOMES (Genome_id), FOREIGN KEY (Genome_id) REFERENCES RAW_ARRAYS (Genome_id) )")
cur.execute("CREATE INDEX IF NOT EXISTS GENOME_ID_INDEX ON SPACER_HITMAP (`Genome_id`)")
cur.execute("CREATE INDEX IF NOT EXISTS Phage_ID_INDEX ON SPACER_HITMAP (`Phage_id`)")
sql_db_connect.commit()


# compute distance scores here
hit_table.to_csv(db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_hitmap.csv",index=False)
# create table phages genomes with 2+ spacer hits.
spacers_2_or_more_master_table = spacers_2_or_more_2.two_or_more(db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_hitmap.csv")
# compute pps distances.
spacer_hitmap_master_table = pps_distance_calculator_annotation_corrected_strand_direction.pps_compute (db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_hitmap.csv" + "_2_or_more_hits.csv")
distance_table = pandas.read_csv(db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_hitmap.csv" + "_2_or_more_hits.csv" + "_distances_annotated.csv")

cur.execute("CREATE TABLE IF NOT EXISTS `PRIMED_SPACER_2+_SPACER_HITMAP_DISTANCES` (Spacer_id STRING,Phage_id STRING,Perc_id STRING,Length STRING, Mismatches STRING,Gapopen STRING,query_start STRING,query_end STRING,Mapped_start_site STRING,Mapped_end_site STRING,evalue STRING,bitscore STRING,Genome_id STRING,orientation STRING,orientation_score STRING,orientation_confidence STRING,questionable_array STRING,array_score STRING,`CRISPR-start` STRING,`CRISPR-end` STRING,repeat_start STRING,repeat_end STRING,spacer_start STRING,spacer_end STRING,dr_repeat_original STRING, dr_repeat_concensous STRING,spacer STRING,Array_tool STRING,array_number STRING,spacer_number STRING,distance STRING,mapped_strand STRING, run STRING, CONSTRAINT Spacer_id PRIMARY KEY (Genome_id, spacer_start, spacer_end), FOREIGN KEY (Genome_id) REFERENCES SPACER_HITMAP(Genome_id), FOREIGN KEY (Phage_id) REFERENCES SPACER_HITMAP(Spacer_GENOME_id)) ")	
cur.execute("CREATE INDEX IF NOT EXISTS GENOME_ID_INDEX ON `PRIMED_SPACER_2+_SPACER_HITMAP_DISTANCES` (Genome_id)")
cur.execute("CREATE INDEX IF NOT EXISTS PHAGE_ID_INDEX ON `PRIMED_SPACER_2+_SPACER_HITMAP_DISTANCES` (Phage_id)")
sql_db_connect.commit()
distance_table.to_sql("PRIMED_SPACER_2+_HITMAP_DISTANCES", sql_db_connect, if_exists='append', index=True)
	
sql_db_connect.close()

# Code to estimate PAMs from mapped sequences.
if (pam_anno == 1):
	filtered_spacer_url = db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered.csv"
	faidx_retriever_genome_block_spacer_mapping_function.mapped_genome_retrieval(filtered_spacer_url,whole_phage_block_url, 20000)
	if not (os.path.isdir(db_directory_path + "phages/")):
		os.mkdir(db_directory_path + "phages/")
	phage_genome_reformatting.reformat(filtered_spacer_url + "_filtered_hits_extracted_faidx_" + "bp_window.fasta", db_directory_path  + "phages/" + b_basename + "_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_genomes.csv") # may want to change the output dir to be in the phages folder
	hits_master_table_url = db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs" + ".lst" + "_reconciled_full_arr_positions.csv", db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_hitmap.csv"
	pam_identification.pam_id(hits_master_table_url, filtered_spacer_url + "_filtered_hits_extracted_faidx_" + "bp_window.fasta",hits_master_table_url + "_pam_summary.csv", 5)
	pam_table = pandas.read_csv(hits_master_table_url + "_pam_summary.csv")
	sql_db_connect = sqlite3.connect("/g/data/va71/crispr_pipeline_annotation/" + "crispr-phage_interactome.sql")
	cur = sql_db_connect.cursor()

	pam_table.to_sql("PAMs",sql_db_connect, if_exists="append",index=False)
	cur.execute("CREATE INDEX IF NOT EXISTS SPACER_ID_INDEX ON PAMS(Genome_id)")
	cur.execute("CREATE INDEX IF NOT EXISTS PHAGE_ID_INDEX ON PAMS(Phage_id)")
	sql_db_connect.close()
if (phage_genomes == 1):
#	code to retrieve mapped sequences without PAM prediction
	if (pam_anno != 1):
		filtered_spacer_url = db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered.csv"
		faidx_retriever_genome_block_spacer_mapping_function.mapped_genome_retrieval(filtered_spacer_url, whole_phage_block_url,20000)
		if not (os.path.isdir(db_directory_path + "phages/")):
			os.mkdir(db_directory_path + "phages/")
		phage_genome_reformatting.reformat(filtered_spacer_url + "_filtered_hits_extracted_faidx_" + "bp_window.fasta", db_directory_path  + "phages/" + b_basename + "_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered.csv") # may want to change the output dir to be in the phages folder


	phage_db_directory_path = db_directory_path + "phages/"
	phage_a =  db_directory_path  + "phages/" + b_basename + "_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered.csv"
	b_phage = filtered_spacer_url + "_filtered_hits_extracted_faidx_" + "bp_window.fasta" 
	phage_master_info_file = open(phage_a + "_summary.csv.txt", "a")
	summary_writer = csv.writer(phage_master_info_file)
	phage_summary_frame = []
	phage_min_seq_id = 0.99
	large_dataset = 0
	phylogeny_switch = 0
	phage_protein_query_switch = 0
	if (phage_protein_query_switch == 1):
		# extract a representative sequence from a query search of genomes/mapped sequences followed by mmseqs2 clustering
		representative_phage_table = genome_representative_clustering.rep_cluster(b_phage, phage_min_seq_id, phage_db_directory_path, phage_a, merge_protein,large_dataset,phylogeny_switch)
		phage_a = representative_phage_table
		b_phage = filtered_spacer_url + "_filtered_hits_extracted_faidx_" + "bp_window.fasta" + "rep_genomes" # need to convert upstream representative clustering to a module
	
	# code to prepare input files for gene annotation
	summary_headers = ["Input_sequence","Input_path","genome_name", "genome_id"]
	phage_summary_frame.append(summary_headers) # This is the phage genome summary table.

	in_seq = b_phage.split("/")
	in_seq = in_seq[-1] # extract just the basename
	genome_sequences = list(SeqIO.parse(b_phage, "fasta")) # these are the phage genomes.

	# convert these genomes into a set!
	phage_dict = {}
	for genome in genome_sequences:
		if (genome.id not in phage_dict):
			phage_dict[genome.id] = genome

	genome_sequences = phage_dict.values()
	sql_db_connect = sqlite3.connect("/g/data/va71/crispr_pipeline_annotation/" + "crispr-phage_interactome.sql")
	cur = sql_db_connect.cursor()
	cur.execute("CREATE TABLE IF NOT EXISTS PHAGE_GENOMES (Phage_id_slice STRING,Phage_id STRING, Phage_id_prefix STRING)")
	k = 0
	while (k < len(genome_sequences)):
		phage_id = genome_sequences[k].id.split("|")
		gene_name = phage_id[0]
		gene_id = phage_id[-1]
		genome_id_prefix = genome_sequences[k].id.split("_")[0]
		res = cur.execute("SELECT Phage_id FROM PHAGE_GENOMES WHERE Phage_id=?",(genome_sequences[k].id,) )
		row = res.fetchone()
		if (row is None):
			cur.execute("INSERT INTO PHAGE_GENOMES (Phage_id_slice,Phage_id, Phage_id_prefix) VALUES (?,?,?)",(genome_sequences[k].id, genome_sequences[k].id.split("-")[0], genome_id_prefix))
			k += 1
		else:
			genome_sequences.pop[k]	

	sql_db_connect.commit()
	genome_sequences = list(SeqIO.parse(b_phage, "fasta")) # these are the phage genomes.
	phage_dict = {}
	for genome in genome_sequences:
		if (genome.id not in phage_dict):
			phage_dict[genome.id] = genome

	genome_sequences = phage_dict.values()
	phage_block_dict = {}
	for gene in genome_sequences:
		gene_id = gene.id.split("|")
		gene_name = gene_id[0]
		gene_id = gene_id[-1]
		phage_block_dict[gene_id] = gene	
		ret_ele = [in_seq, phage_a, gene_name, gene_id]
		phage_summary_frame.append(ret_ele)
	
	# construct genome summary csv file for mapped sequences
	for frame in phage_summary_frame:
		summary_writer.writerow(frame) 
	phage_master_info_file.close()
	SeqIO.write(genome_sequences,b_phage, "fasta") 
	
	# concatenate to a file containing an indexed form of all mapped phage sequences
	if (os.path.isfile(phage_hits_block + ".fai")):
		# need to only add genomes which are not already in the data block. Extract the ids from the samtools file
		with open (phage_hits_block + ".fai", "r") as csvfile:
			fai_table = list(csv.reader(csvfile,delimiter='\t'))
		subprocess.run(["rm " + phage_hits_block + ".fai"],shell=True)	
		fai_headers = set()
		for row in fai_table:
			if (row[0] not in fai_headers):
				fai_headers.add(row[0])
			else: 
				print("Error!! Should not have duplicate fai headers in one job!!")
		new_genomes = SeqIO.parse(genome_sequences, "fasta")
		existing_block = open(phage_hits_block, "a")
		for genome in genome_sequences:
			if (genome.id not in fai_headers):
				SeqIO.write(genome, existing_block, "fasta")
		existing_block.close()
	else:	
		#	print(b)
		subprocess.run(["cat " + b_phage + " >> " + phage_hits_block],shell=True)
	
	subprocess.run(["/g/data/va71/crispr_pipeline_annotation/seqkit faidx " + "--update-faidx " + phage_hits_block],shell=True)

if (protein_generation_switch == 1):
	input_genome_file_name = b_phage
	generate_protein(input_genome_file_name,phage_protein_hits_block, rnafold_switch, merge_protein)

# start of gene_annotation_summary_table_generation for host encoded sequences

with open (a, "r") as csvfile:
	hit_table = list(csv.reader(csvfile))
sequence_master_info_file = open(a + "_summary.csv.txt", "w") # or could just use print statements!! Nah, may generate multiple files. Worth opening an independent stream.
# need to create a table specifying the protein query, matching genome, extracted all orf positions and identity of each orf using HHpred/HMMsearch
summary_writer = csv.writer(sequence_master_info_file)
summary_frame = []
if (protein_query_switch == 1):
	b = b + "_rep_genomes.fasta" #switch to representative genomes of a clustered representation of the host set.
# may want to make this summary file a option in the CLI!! Often not needed for pipeline running
# summary table file generation for when an input query sequence was used as input
if (db_directory_switch != 0 ): # If a csv/query sequence is used for input
	summary_headers = ["Input_sequence","Input_table_name","Input_table_entries","genome_file_names","corresponding_genome_ids", "genome_protein_orfs"]
	summary_frame.append(summary_headers)
# variable names:
	in_seq = b.split("/")
	in_seq = in_seq[-1] # extract just the basename
	for row in hit_table:
		ret_ele = [in_seq, a, row[3]] # 
		summary_frame.append(ret_ele)
	genome_sequences = list(SeqIO.parse(b, "fasta"))
	block_dict = {}
	for gene in genome_sequences:
		gene_id = gene.id.split("|")
		gene_id = gene_id[-1]
		block_dict[gene_id] = gene
	i = 1
	while (i < len(summary_frame)):
		genome_row = summary_frame[i][2].split("|") 
		gene_id = str(genome_row[1])
		gene_id_suffix = gene_id.split("::")
		gene_id_suffix2 = gene_id_suffix[1]
		gene_id_suffix3 = gene_id_suffix2.split(":")
		gene_id_suffix4 = gene_id_suffix3[1]
		gene_id_suffix5 = gene_id_suffix4.split("_")
		gene_id_suffix5 = gene_id_suffix5[0]
		true_gene_id = gene_id_suffix[0] + "::" + gene_id_suffix3[0] + ":" + gene_id_suffix5
		if (true_gene_id) in block_dict: # if a given genome_id is in the table
			genome = block_dict[true_gene_id]
			if (protein_query_switch == 1):
				genome_id = genome.id  
				genome_file_name = genome.description.split("|")
				genome_file_name = genome_file_name
			else:
				genome_id = genome.description.split("|")
				genome_file_name = genome_id[0]
				genome_file_name = genome_file_name.split(" ")
				genome_file_name = genome_file_name[0] # need to switch the genome file and genome id variable names!!
				genome_id = genome_id[1]
			summary_frame[i].extend([genome_id, genome_file_name, b + "_output_full.txt"]) # protein orf_name
		else:
			print("Error!!")
		i += 1
	for frame in summary_frame:
		summary_writer.writerow(frame) 
	sequence_master_info_file.close()
else: #  case for when host genomes only are used as input
	summary_headers = ["Input_sequence","Input_path","genome_name", "genome_ids"]
	summary_frame.append(summary_headers)
	in_seq = b.split("/")
	in_seq = in_seq[-1] # extract just the basename
	genome_sequences = list(SeqIO.parse(b, "fasta")) # not sure how well this command will work at scale!!!
	block_dict = {}
	for gene in genome_sequences:
		gene_id = gene.id.split("|")
		gene_name = gene_id[0]
		gene_id = gene_id[-1]
		block_dict[gene_id] = gene	
		ret_ele = [in_seq, a, gene_name, gene_id]
		summary_frame.append(ret_ele)
	for frame in summary_frame:
		summary_writer.writerow(frame) # could do this concurrently with writing in memory!!
	sequence_master_info_file.close()

if (phage_genomes == 1):
	my_pickle_phage = open(phage_a + "_pickle_block_dict.pickle", "wb")
	pickle.dump(phage_block_dict, my_pickle_phage)
	my_pickle_phage.close()
else:		
	my_pickle = open(a + "_pickle_block_dict.pickle", "wb")
	pickle.dump(block_dict, my_pickle)
	my_pickle.close()

protein_annotation_frames = []

# bypass switch in the case of spacers already being mapped!!
if (bypass_switch == 1):
	b_basename = b.split("/")
	b_basename = b_basename[-1]
	cores = 48
	# block directory
	perc_identity = 0.85
	input_spacer_fasta = db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs" + ".lst" + "_full_real_arr_positions.csv" + "_spacers.fasta"	
	spacer_dict = spacer_dict_generator.spacer_dict_gen(input_spacer_fasta + "_all_hits.csv")
	

# main function calls for gene annotation in parallel:
# Note this was run as a seperate standalone workflow, partially because the compute requirements for both workflows were too large for a single batch run (for larger subtypes)
if (phage_genomes == 1):
	phage_protein_annotation_frame = []
	phage_summary_frame = phage_summary_frame[1:] # may not need this if phage summary frame already lacks headers
	pool = Pool(cores)
	db_directory_switch = 0
	target_pos_switch = 0
	for frame in phage_summary_frame:
		pool.apply_async(annotation_parallelisation,(frame, phage_block_dict, db_directory_switch, b_phage, phage_db_directory_path, phmmer_switch,hhsearch_switch,hhblits_switch,dali_switch,target_pos_switch,pfam_switch, criscasfinder_switch, merge_protein,eggnog_switch,vmatch_switch,virsorter_switch,rnafold_switch,padlocplus_switch_novel, True))
	pool.close()
	pool.join()


print("checkpoint 5")
print(db_directory_switch)
print(db_directory_path)
print(b)
summary_frame = summary_frame[1:]
pool = Pool(cores)
for frame in summary_frame:
	pool.apply_async(annotation_parallelisation, (frame, block_dict, db_directory_switch, b, db_directory_path, phmmer_switch,hhsearch_switch,hhblits_switch, dali_switch, target_pos_switch, pfam_switch, criscasfinder_switch, merge_protein,eggnog_switch, vmatch_switch,virsorter_switch,rnafold_switch,padlocplus_switch_novel, False))
pool.close()
pool.join()
print("Done!")
