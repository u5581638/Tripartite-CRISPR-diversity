

import csv
import sys
import sqlite3
from Bio import SeqIO
import os
import subprocess

# script to run PlasME for each genome/ mapped sequence in a standalone manner
# INPUT: 1. path to folder containing DNA extracted within 20kb of CRISPRs (db_directory_path)
#		 i.e. /g/data/va71/crispr_pipeline_annotation/annotation_upgraded_main_workflow_run/cas13a/genomes/
#		 2. query search name 
#		 i.e. cas13a.fasta
#		 3. folder containing mapped target files (phage_db_directory_path)
#	 	 i.e. /g/data/va71/Alex/crispr_pipeline_annotation/annotation_upgraded_main_workflow_run/cas13b/genomes/phages/
#		 4. mapped target contigs from spacer mapping (in FASTA format)
#		 i.e. cas13b.fasta_all_hits.csv_genomes.fasta_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_hitmap.csv_standardised.csv_non_redundant.csv_filtered.csv_deduplicated.csv_filtered_hits_extracted_faidx_bp_window.fasta
#		 5. number of CPUs to use 
#		 i.e. 24
# OUTPUT: csv file containing PLASMe scores for linear/circular DNA.
# SHELL: see ...

phage_genomes = 1
#db_directory_path = "/g/data/va71/crispr_pipeline_annotation/annotation_upgrades_test_sequences/cas12a_test2_4/genomes/"
db_directory_path = sys.argv[1]
#b_basename = "cas12_archtype.fasta"
b_basename = sys.argv[2]
filtered_spacer_url = db_directory_path + "spacer_distribution_analysis/" + b_basename + "_all_hits.csv_genomes.fasta" + "_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_hitmap.csv" + "_standardised.csv_non_redundant.csv_filtered.csv"
phage_db_directory_path = db_directory_path + "phages/"
phage_a =  db_directory_path  + "phages/" + b_basename + "_all_hits.csv_genomes.fasta" + "_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered.csv"
b_phage = sys.argv[3] 

# cores = 24
cores = sys.argv[4]
# a = "/g/data/va71/crispr_pipeline_annotation/annotation_upgrades_test_sequences/cas12a_test2_4/queries/cas12_archtype.fasta_all_hits.csv"
#a = sys.argv[5]
#b = a + "_genomes.fasta"

subprocess.run(["/g/data/va71/crispr_pipeline_annotation/plasme.sh " + b + " " + b + "_p_ovrlps.fa " + str(cores) + " " + db_directory_path],shell=True)
subprocess.run(["/g/data/va71/crispr_pipeline_annotation/plasme.sh " + b_phage + " " + b_phage + "_p.fa " + str(cores) + " " + phage_db_directory_path],shell=True)