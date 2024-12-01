
# protein annotation workflow api
import sys
if ('/g/data/va71/crispr_pipeline_annotation' not in sys.path):
	sys.path.append('/g/data/va71/crispr_pipeline_annotation')

import sqlite3
from Bio import SeqIO
import csv
import os
import subprocess
import pilercr_pos_extractor_annotation_full_spacers
import shutil
import timeit
from time import process_time
from time import time
import fnmatch
import pyfasta_but_faster
from multiprocessing import Pool
import pickle
# from annotation_parallelisation_function_SQL8 import annotation_parallelisation

import genemark_sequence_reformatting
import prodigal_genemark_reconciliation
import genome_representative_clustering
from datetime import date
import pandas
import random
import spacer_mapping_parallelisation

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
#	print(blast_db_seqs)
	print("loop start!!")
	i=0
	genome_path = blast_db_seqs[0]
	for all_seqs in blast_db_seqs:
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
			
	return genome_path

# helper function to run makeblastdb for db generation
def db_generation_subroutine (db_directory_path):
	# need to format the db fasta file sequences with their respective block names
#	print("db_file_is: " + db_directory_path)

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
		subprocess.run(["samtools", "faidx", input_genome_file_name + "_CDS_raw.fnn"])
	else:		
		subprocess.run(["/g/data/va71/crispr_pipeline_annotation/prodigal/prodigal", "-p", "meta", "-i", input_genome_file_name, "-o", input_genome_file_name + "_output_full.txt", "-a", input_genome_file_name + "_aa_raw.fasta" ]) # output should be redirected to genome specific folder
		if (merge_protein == 1):
				subprocess.run(["/g/data/va71/crispr_pipeline_annotation/annotation_database/annotation_upgrades/protein_reconciliation/genemark_reinstall/gms2_linux_64/gms2.pl", "--seq", input_genome_file_name,"--genome-type", "auto", "--out", input_genome_file_name + "_trans_aa_gms2.txt", "--faa", input_genome_file_name + "_trans_aa_gms2.faa", "--max-iter", "3","--conv-thresh", "0.9"])
				genemark_sequence_reformatting.genemarkS2_reformatting(input_genome_file_name + "_trans_aa_gms2.faa")
				prodigal_genemark_reconciliation.union(input_genome_file_name + "_trans_aa_gms2.faa" + "_geneS2prod_reformatted.fasta", input_genome_file_name + "_aa_raw.fasta", input_genome_file_name + "_aa_raw.fasta") # This will hopefully override the original protein file
	# This is where samtools is called. Need to add a function on input_protein_file to replace the whitespace with an allowed character.
#	hash_adder(input_genome_file_name + "_aa_raw.fasta") # replace the whitespace with #
	subprocess.run(["/g/data/va71/crispr_pipeline_annotation/seqkit", "faidx", input_genome_file_name + "_aa_raw.fasta"])	
	subprocess.run(["makeblastdb", "-in", input_genome_file_name + "_aa_raw.fasta", "-dbtype", "prot"])
	subprocess.run(["cat " + input_genome_file_name + "_aa_raw.fasta >> " + protein_block_name],shell=True)
	# may need to improve this to
	# rmdup -i " + genome_hits_block + " -o " + genome_hits_block  + " -d duplicated.fa.gz -D duplicated.detail.txt"
	subprocess.run(["/g/data/va71/crispr_pipeline_annotation/seqkit " + "rmdup -i " + protein_block_name + " -o " + protein_block_name  + " -d duplicated.fa.gz -D duplicated.detail.txt"],shell=True)
	subprocess.run(["samtools", "faidx", protein_block_name])
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

# generate table compatible with SQL. Note: this function was not used.
def genome_hit_sql_table_generation(hit_table, genome_url, params, sql_db_connect):
#	con = sqlite3.connect(sql_db_connect)
#	cur = con.cursor()
	ret_list = []
	header = ["RUN","Seed_id", "Genome_id", "perc_identity", "length", "mismatch", "gapopen", "query_start", "query_end", "match_start", "match_end","evalue","bitscore", "my_date", "params","Genome_id_prefix","Genome_path"]
	ret_list.append(header)
	out_ret = open(genome_url, "w")
	spam_writer = csv.writer(out_ret)
	spam_writer.writerow(header)
	for hit in hit_table:
		genome_id_prefix = hit[1].split("_") [0]
		if (genome_id_prefix[:1] == "Ga"):
			row = ["Not_assigned"] + hit + [date.today(), params] + [genome_id_prefix]
		else:
			row = ["Not_assigned"] + hit + [date.today(), params] + [None] # These may need to be null values	
		spam_writer.writerow(row)
		ret_list.append(row)
	out_ret.close()	
	return 0	

# Is there a simpler way of making a CLI? See the run_alphafold code!!
# Retrieve genome sequences by identifier.
def genome_hit_lookup (hit_table, db_directory_path, genome_file_url, block, genome_row_no=2): # b?
	row_dict = {}
	if (os.path.isfile(genome_file_url)):
		os.remove(genome_file_url)
	genome_file = open(genome_file_url, "a")
	# print(hit_table)
	for row in hit_table:
		# go through hit table and split based on the genome identifier
		genome_row = row[genome_row_no]
		genome_row = genome_row.split("|") 
		true_gene_id = genome_row[0]
		'''
		gene_id = str(genome_row[1])
		gene_id_suffix = gene_id.split("::")
		gene_id_suffix2 = gene_id_suffix[1]
		gene_id_suffix3 = gene_id_suffix2.split(":")
		gene_id_suffix4 = gene_id_suffix3[1]
		gene_id_suffix5 = gene_id_suffix4.split("_")
		gene_id_suffix5 = gene_id_suffix5[0]
		true_gene_id = gene_id_suffix[0] + "::" + gene_id_suffix3[0] + ":" + gene_id_suffix5		
		'''
		# only need to construct a new list/dict if the current genome_row genome differs from the previous!	
		genome_row = block # in theory this should be the block
		if (genome_row not in row_dict):
			my_id = {true_gene_id}
			row_dict[genome_row] = my_id
		else:
			row_dict[genome_row].add(true_gene_id) 
	for row in row_dict.items():
		genome_row = row[0]
		block = SeqIO.parse(genome_row, "fasta")
		for gene in block:
			gene_key = gene.id.split("|")
		#	true_gene_id = 
			if (len(gene_key) == 1):
			#	print("mate!!")
				gene_key = gene_key[0]
			else:
				# this only works when the input data is either window or global genome blocks - not generalisable!!
				gene_key = gene_key[1]	
			if (gene_key in row[1]):
				my_genome = gene
				my_genome.id = gene_key
				my_genome.description = ""  # my_genome.description # + "|" + row[0] # now sequence should be annotated with block + genome_id + sequence_id
				# need to modify this code to erase the block name!! This is probably the easiest fix!!
				SeqIO.write(my_genome, genome_file, "fasta") # This needs to be written to file then indexed with samtools
	genome_file.close()
	return 0

#START!!
# code for directly running from cmdline goes here!!
# This is the main code for performing CRISPR array validation, which is followed by spacer mapping and gene annotation
# This was originally designed to run as a CLI program with diverse functions. However, this was not implemented and a single workable workflow was run instead.

# input_flags:
# -i input (either csv_file or input sequence in fasta format)
# -t input_type (query/table)
# -p use phmmer
# -b use hhblits
# -s use hhsearch
# -a use alphafold + DALI
# should supply a number of potential tokens as well as having default values. Tokens are represented with an '-' specifier.

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
skip_genome_dup_detection = 0 



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
		case '-skip_genome_dup_detection':
			skip_genome_dup_detection = 1
			i += 1
						

	print(i)	

eggnog_switch=0
vmatch_switch=0
virsorter_switch=0
rnafold_switch=0
padlocplus_switch_novel=0			


print("checkpoint 0")

# sql_db_connect = sqlite3.connect("/g/data/va71/crispr_pipeline_annotation/" + "crispr-phage_interactome.sql")
sql_db_connect = sqlite3.connect("/g/data/va71/crispr_pipeline_annotation/" + "crispr-phage_interactome_main_run2.sql")
genome_hits_block = "/g/data/va71/crispr_pipeline_annotation/" + "genome_annotation_block_main_run.fasta" 
phage_hits_block = "/g/data/va71/crispr_pipeline_annotation/" + "phage_annotation_main_run.fasta"
protein_hits_block = "/g/data/va71/crispr_pipeline_annotation/" + "protein_annotation_main_run.fasta"
phage_protein_hits_block = "/g/data/va71/crispr_pipeline_annotation/" + "phage_protein_annotation_main_run.fasta"
whole_phage_block_url = "/g/data/va71/crispr_pipeline_annotation/" + "whole_phage_genome_block_main_run.fasta"
cur = sql_db_connect.cursor()

golddata_url = "/g/data/va71/crispr_pipeline_annotation/golddata/golddata_7_6_2023.xlsx"
if (spacer_generation_bypass_switch == 0):
	if (gold_switch == 1):
		gold_study = pandas.read_excel(golddata_url,sheet_name='Study')
		gold_biosample = pandas.read_excel(golddata_url,sheet_name='Biosample')
		gold_Organism = pandas.read_excel(golddata_url,sheet_name='Organism')
		gold_Sequencing_Project = pandas.read_excel(golddata_url,sheet_name='Sequencing Project')
		gold_Analysis_project = pandas.read_excel(golddata_url,sheet_name="Analysis Project")


	# Here I can get away with zero constraints because NCBI already ensures table consistency. In all other sql tables though constraints are required. Will need to add these contraints in a create table statements
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


	# now should either an input csv table, or csv table generated from the query sequence.

	# also want to create a version of the protein capable of non-specific annotation given an input set of genomes (db_directory_switch_option == 0)
	print("checkpoint 1!!")

	# need to retrieve the corresponding labelled genome.
	# Is this neccessary?
	print("db_directory_switch:")
	print(db_directory_switch)
	if (db_directory_switch != 0):
		if (input_type == "query"): # here we want to take the table generated from db_generation
			mod_file_path = file_path # .split('/')
		#	mod_file_path = mod_file_path[-1]
			a = mod_file_path + "_all_hits.csv" #  This should point to the generated table. May need a '/' between db_directory_path and file path for this to work!!
		elif (input_type == "table"): # here we want to take the input table, which has been labelled!!
			a = mod_file_path # this will be a specified table.
		else: 
			print("problem with the input character!!")



		with open (a, "r") as csvfile:
			hit_table = list(csv.reader(csvfile))
			#print(hit_table)
		# for each sequence in the query file, find the corresponding genome hit
		# need to think about the best way to do this!!
		# probably want to annotate with 
		b = a + "_genomes.fasta" # need to open genomes!!!!! this is the real problem
		genome_hit_lookup(hit_table, db_directory_path, b, genome_path, 1)
		# bug must be past this block!!
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
		#	print(fai_table)
			for row in fai_table:
				if (row[0] not in fai_headers):
					fai_headers.add(row[0])
				else: 
					print("Error!! Should not have duplicate fai headers in one job!!")
			new_genomes = SeqIO.parse(genome_file, "fasta")
			existing_block = open(genome_hits_block, "a")
			for genome in new_genomes:
				if (genome.id not in fai_headers):
					SeqIO.write(genome, existing_block, "fasta")
			existing_block.close()


		else:	
		#	print(b)
			subprocess.run(["cat " + b + " >> " + genome_hits_block],shell=True)
		subprocess.run(["/g/data/va71/crispr_pipeline_annotation/seqkit rmdup -i " + genome_hits_block + " -o " + genome_hits_block  + " -d duplicated.fa.gz -D duplicated.detail.txt"],shell=True)
		print("Passed stdout bottleneck!!")
		subprocess.run(["/g/data/va71/crispr_pipeline_annotation/seqkit faidx " + "--update-faidx " + "\'" + genome_hits_block + "\'"],shell=True)
		# if option -q is used then prodigal should be used to predict proteins arising from tblastn hits.
		# blastp should then identify these proteins.
		# these proteins should then be retrieved!!
		# these should also be a set, so that redundant proteins are removed.
		if (protein_query_switch == 1):
			representative_table = genome_representative_clustering.rep_cluster(b, min_seq_id, db_directory_path, a, merge_protein, large_dataset, phylogeny_switch, sys.argv)
			a = representative_table
			with open (a, "r") as csvfile:
				hit_table = list(csv.reader(csvfile))

		genome_hit_sql_table_generation(hit_table,a, "".join(sys.argv), sql_db_connect) # need to reach here to bypass creating new genome blocks. The reason is because this function modifies the hit_table!!!
		print("Good up to here?")
		# RUN this code as a standalone. Attempt to integrate the prefix ids. Should not need to rest!!!!
		export_table = pandas.read_csv(a,header=0)
		genome_prefix_ids = {}
		genome_type = {}
		genome_is = "HOST"
	#	print(export_table.columns)
		for entry in export_table['Genome_id']:
			my_entry = entry
			my_entry = str(my_entry).split("_") 	
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
		

		export_table = export_table.assign(Genome_id_prefix = list(genome_prefix_ids.values()) [0]) # These don't seem to work!
		export_table = export_table.assign(Genome_type = list(genome_type.values()) [0]) # These don't seem to work!
		export_table = export_table.assign(Genome_path = genome_path)

	#	cur.execute("PRAGMA foreign_keys = ON;")
	#	cur.execute('CREATE TABLE IF NOT EXISTS GENOMES (RUN STRING, Seed_id STRING, Genome_id STRING, perc_identity STRING, length STRING, mismatch STRING, gapopen STRING, query_start STRING, query_end STRING, match_start STRING, match_end STRING, evalue STRING, bitscore STRING, my_date STRING, params STRING, Genome_type STRING, Genome_id_prefix TEXT, Genome_path TEXT)') # FOREIGN KEY (Genome_id_prefix) REFERENCES GOLD_ANALYSIS_PROJECT (`AP GOLD ID`) see if this works without the foreign key. Might not be feasible to make mapping 1 to 1 due to NCBI IDs.
		sql_db_connect.commit()
		# generate a unique key or keep trying
		while (True):
			index = random.randint(1000000000000000000, 10000000000000000000)
			index = str(index)
			break
		#	index = random.choices(L, k = 1) [0]
		#	res = cur.execute('SELECT * FROM GENOMES WHERE RUN LIKE ?', (index,))
		#	row = res.fetchone()
		#	if (row != None): # There should not be an identical run number. THIS MAY BE A SOURCE OF BUGS!!!!!
		#		print("random_generation!!") # This should
		#	else:
		#		break
					
		i = 0
		# should make a standalone code file containing just this loop and the nessessary inputs.
		# Try just one instance of export table so that each command can be run interactively
		# It might be possible to make an optimisation here because, in theory, each genomes table only contains one seed!!!!!!!!
		# This would reduce the time complexity to a single O(n) pass.
		'''
		if (skip_genome_dup_detection == 0):
			while i < len( export_table):
				res = cur.execute('SELECT * FROM GENOMES WHERE Genome_id LIKE ?', ('%'+ export_table.iloc[i][2] + '%',)) # hang on, do these need to be genomes????? Matching between Genomes and seeds = wrong!!!!!	
				row = res.fetchone()
				if (row is not None):
					# In this instance need to add seed sequence to first row.  (.
				#	print(row[0])
				#	print(export_table.iloc[i][1])
					new_value = str(row[1]) + "," + export_table.iloc[i][1]
					cur.execute('UPDATE GENOMES SET Seed_id =? WHERE Genome_id LIKE ?', (new_value, '%' + export_table.iloc[i][2] + '%',))			# need to create a table if it does not already exist mapping genome ids to seed values
				#	print(export_table)
					export_table.drop(labels=export_table.index[i],axis='index',inplace=True) # Remove overlapping genome from the collection of genomes to annotate. Might still want to spacer map though???? May need to create a seperate set of indexed files and tables to handle this!! 
				#	print(export_table)
				#	export_table.iloc[0][i] = index
				else:
					# This is the problematic line!!
					export_table.iat[i, 0] = index
					i += 1
				
				sql_db_connect.commit()	
		'''
			# Need to decide whether to keep this dumpster-fire of sql code!!
			#	cur.execute('SELECT * FROM GENOME_SEED_JUNCTION')
			#	names = list(map(lambda x: x[0], cur.description))
			#	previous = names[-1].split("_") [-1]
			#	new = previous += 1
			#	cur.execute('ALTER TABLE GENOME_SEED_JUNCTION ADD ? STRING', ("seed_" + str(new))) # not sure if this operation is allowed in sqlite3
			#	cur.execute('UPDATE GENOME_SEED_JUNCTION SET seed_?=? WHERE Genome_id=?;', (str(new),export_table[i][0],export_table[i][1])) # need to check with someone who knows SQL whether this makes sense/could be written better!! Fortunately this code is sort of optional!!
				
			
				# Answer: Spacer mapping does not need to occur again PROVIDED previously mapped phage genomes are associated with the new genomes. MAKE SURE THIS IS DONE BEFORE RUNNING THE PIPELINE!!!!
				
		#	else:
				# passively add to seed table. A lot more complicated then changing the field
			#	cur.execute('CREATE TABLE IF NOT EXISTS GENOME_SEED_JUNCTION (Genome_id STRING PRIMARY KEY)')
			#	cur.commit()


			#	cur.execute('SELECT * FROM GENOME_SEED_JUNCTION')
			#	names = list(map(lambda x: x[0], cur.description))
			#	previous = names[-1].split("_") [-1]
			#	new = previous += 1
			#	if (previous == "id"):
			#		cur.execute('ALTER TABLE GENOME_SEED_JUNCTION ADD ? STRING', ("seed_0"))
			#		cur.execute('INSERT INTO GENOME_SEED_JUNCTION(Genome_id,seed_0) VALUES(?,?)', (export_table[i][1],export_table[i][0]))
			#	else:
			#		cur.execute('ALTER TABLE GENOME_SEED_JUNCTION ADD ? STRING', ("seed_" + str(new)))
			#		cur.execute('INSERT INTO GENOME_SEED_JUNCTION(Genome_id,seed_? VALUES(?,?)', (str(new),export_table[i][1],export_table[i][0])) # need to check with someone who knows SQL whether this makes sense/could be written better!! Fortunately this code is sort of optional!!

			#	cur.commit()
			# end of sql dumpster-fire code!!
	#	cur.execute('CREATE INDEX IF NOT EXISTS GENOME_PREFIX_INDEX ON GENOMES(Genome_id_prefix)')
	#	sql_db_connect.commit()
		# failure adding genomes
	#	export_table.to_sql('GENOMES',sql_db_connect, if_exists='append', index=False)		
		# NEED TO DECIDE WHETHER TO INCLUDE THE LINE BELOW> THE ABSENCE OF HITS IS BECAUSE EACH ENTRY ALREADY EXISTS IN THE SQL TABLE!!	
	#	export_table.to_csv(a)
		# 

		# need to sort hit table such that hits from the same database file are searched in a single pass of the database
		# maybe instead step 1 should be to dictionalise the rows with the same common block (key)
		# for each entry:
		# load a different block
		# if block id and row id matches then write to genome_file
		# continue iterating through all the blocks	
		# genomes should be accessible by index!!
	#	print ("Done IO!!")
		

	#	print(genome_file)

	# this step is all about extracting genomes

	else: # need an if statement to delete this directory prior to running the program. Ditto for BLAST compiled db!!
		# Save fixing this for last!! Don't need to run under default circumstances!!!!!!!!!!
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
			#	genome_prefix_ids['Genome_id_prefix'] = [my_entry]
				genome_type['Genome_type'] = [genome_is]
			else:
			#	genome_prefix_ids.append(my_entry)
				genome_type.append(genome_type)
			
	#	export_table['Genome_id_prefix'] = genome_prefix_ids.values()
		export_table['Genome_type'] = list(genome_type.values()) [0]
	#	i = 0
	#	while i < len( export_table):
	#		res = cur.execute('SELECT * FROM GENOMES WHERE Genome_id=?', tuple(hit[1]))
			
	#		if (res.fetchone() is not None):
				# In this instance just remove entry from export table.
	#			export_table.drop(labels=i,axis=0,inplace=True)
	#		else:
	#			i += 1	

		cur.execute('CREATE TABLE IF NOT EXISTS GENOMES (RUN STRING, Seed_id STRING, Genome_id STRING, perc_identity STRING, length STRING, mismatch STRING, gapopen STRING, query_start STRING, query_end STRING, match_start STRING, match_end STRING, evalue STRING, bitscore STRING, my_date STRING, params STRING, Genome_type STRING, Genome_id_prefix STRING, Genome_path TEXT, FOREIGN KEY (Genome_id_prefix) REFERENCES GOLD_ANALYSIS_PROJECT (`AP GOLD ID`))')
		cur.execute('CREATE INDEX IF NOT EXISTS GENOME_PREFIX_INDEX ON GENOMES(Genome_id_prefix)')
		cur.commit()
		export_table.to_sql('GENOMES',sql_db_connect, if_exists='append', index=True)


		SeqIO.write(my_genomes, b, "fasta") # rename the genome file to be consistent with the other settings. Better to use rename?

		# need to create an input table directly from the genome files encoding the same information!!
		# need to standardise the table input for both cases in sqlGe
	index_file = open(file_path + "_index.txt", "w")
	index_file.write(str(index))
	index_file.close()
	# need to add primary and foreign keys to sql table!!
	sql_db_connect.close()
	# now need to consider garbage collection

	# mistake with head genome retrieval must be before this line!!!!!!!!!!!!!!!!!!!!!!!!
	print("checkpoint 2!!")

	# write prodigal + genemark orf translation of target genomes here!!
	# index with samtools!!
	if (protein_generation_switch == 1):
		input_genome_file_name = b
		generate_protein(input_genome_file_name,protein_hits_block, rnafold_switch, merge_protein)
	print("checkpoint 3!!")
	print(db_directory_path)
	#cores = 48  this should be a command line option

	if (spacer_mapping_switch == 1):
		print("checkpoint 4!!")
		if (os.path.isdir(db_directory_path + "spacer_distribution_analysis/")):
				pass 
		else:	
			os.mkdir(db_directory_path + "spacer_distribution_analysis/")
		b_basename = b.split("/")
		b_basename = b_basename[-1]
		print (b)
		print(b_basename)
		# need to seperate split spacers into smaller files
		output_dir = db_directory_path + "spacer_distribution_analysis/" + b_basename
		partition_size = os.path.getsize(b) // cores
		partitions = pyfasta_but_faster.rename(b,output_dir,  1966080) # was 1966080
		m = 0
		# may be able to parallelise this loop to save KSU!!
		partition_cores = 1
		protopool = Pool(int(cores // partition_cores))
	
	# function to run CRISPR-array prediction and validation.	
	# This is parallelised to compensate for the dramatically increased run time due to soo many prediction tools running.
		while (m <= partitions):
			partition_dir = output_dir + "_partition_" + str(m) + ".fasta"
			protopool.apply_async(spacer_mapping_parallelisation.parallel_spacer_partition, (output_dir,partition_dir, db_directory_path,b_basename,crispr_detect, crispr_orientation, partition_cores))
			m += 1
		protopool.close()
		protopool.join()
else:
	b_basename = file_path + "_all_hits.csv" + "_genomes.fasta"
	b_basename = b_basename.split("/") [-1]
	# Need to put in a bypass to jump to this line of code!!!!!!

# table_reverse(db_directory_path + "spacer_distribution_analysis/" + b_basename + "_crisprs" + ".lst" + "_reconciled_full_arr_positions.csv")

# Merge CRISPR-array predictions
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

	# now need to blast against an input db
	# may need to define parameters as optional flags!!
	# need to launch script with maximum cpus.
	# block directory	
exit()
	# this is end of spacer_preproccessing!!

