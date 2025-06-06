import csv
import sys
import sqlite3
from Bio import SeqIO
import os
import subprocess
from multiprocessing import Pool
import pickle
from annotation_parallelisation_function_SQL11 import annotation_parallelisation
import virsorter_parser2

# got to get the phage block dict!!
# get phage_block_dict_from_pickle!!
# need to retrieve/load phage_summary_frame!!
# also need to check exactly how proteins from generate protein are handled!!!

# script to annotate contigs from each subtype:
# INPUT: 1. path to folder containing DNA extracted within 20kb of CRISPRs (db_directory_path)
#		 i.e. /g/data/va71/crispr_pipeline_annotation/annotation_upgraded_main_workflow_run/cas13b/genomes/
#		 2. query search name 
#		 i.e. cas13b.fasta
#		 3. tBLASTn result table from contig retrieval.
#		 i.e. /g/data/va71/crispr_pipeline_annotation/annotation_upgraded_main_workflow_run/cas13b/queries/cas13b.fasta_all_hits.csv
#		 4. table containing mapped spacer hits (filtered_spacer_url)
#		 i.e. cas13b.fasta_all_hits.csv_genomes.fasta_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_hitmap.csv_standardised.csv_non_redundant.csv_filtered.csv
#        5. folder containing mapped target files (phage_db_directory_path)
#	 	 i.e. /g/data/va71/Alex/crispr_pipeline_annotation/annotation_upgraded_main_workflow_run/cas13b/genomes/phages/
#		 6. summary file of mapped target contigs (phage_a)
#		 i.e. cas13b.fasta_all_hits.csv_genomes.fasta_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered.csv_summary.csv.txt
#		 7. mapped target contigs from spacer mapping (in FASTA format)
#		 i.e. cas13b.fasta_all_hits.csv_genomes.fasta_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_hitmap.csv_standardised.csv_non_redundant.csv_filtered.csv_deduplicated.csv_filtered_hits_extracted_faidx_bp_window.fasta
#		 8. Additionally, this script depends on the output and directories used to perform spacer mapping (see "annotation_cas13a_working_spacer_mapping1.sh" as an example) 
# OUTPUT: 1. annotations of each gene in each contig using Pfam, PDB70 and DEFLOC via HMM searches.
#		  2. Virsorter predictions (depending on the exact command used)
#		  3. Viennafold RNA secondary structure predictions of each ORF employed. (depending on the exact parameters/switches used in input section)

# SHELL: see cas13b_annotation_only_plasme_general_run.sh or cas13b_annotation_phage_only_plasme_general_run.sh or cas13b_virsort_only_plasme_general_run_all.sh for different run parameters
phage_genomes = 1
# db_directory_path = "/g/data/va71/crispr_pipeline_annotation/annotation_upgrades_test_sequences/cas12a_test2_4/genomes/"
db_directory_path = sys.argv[1]

#b_basename = "cas12_archtype.fasta"
b_basename = sys.argv[2]
filtered_spacer_url = db_directory_path + "spacer_distribution_analysis/" + b_basename + "_all_hits.csv_genomes.fasta" + "_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_hitmap.csv" + "_standardised.csv_non_redundant.csv_filtered.csv"
phage_db_directory_path = db_directory_path + "phages/"
phage_a =  db_directory_path  + "phages/" + b_basename + "_all_hits.csv_genomes.fasta" + "_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered.csv"
b_phage = filtered_spacer_url + "_deduplicated.csv" + "_filtered_hits_extracted_faidx_" + "bp_window.fasta" 
db_directory_switch = 0
phmmer_switch = 0 
hhsearch_switch = 0
hhblits_switch = 1
dali_switch = 0
target_pos_switch = 0
pfam_switch = 1
criscasfinder_switch = 1
merge_protein = 1
cores = int(sys.argv[4])
vmatch_switch = 0
virsorter_switch = int(sys.argv[5])
rnafold_switch = 1
padlocplus_switch_novel = 0
eggnog_switch = 1
phage_only = 0
host_only = 1

if (virsorter_switch != 1):
	phage_master_info_file = open(phage_a + "_summary.csv.txt", "r")
	phage_summary_frame = list(csv.reader(phage_master_info_file))
	phage_master_info_file.close()


#a = "/g/data/va71/crispr_pipeline_annotation/annotation_upgrades_test_sequences/cas12a_test2_4/queries/cas12_archtype.fasta_all_hits.csv"
a = sys.argv[3]
b = a + "_genomes.fasta"



if (virsorter_switch == 1):
	subprocess.run(["source /g/data/va71/my_conda/new_conda/etc/profile.d/conda.sh && conda activate vs2 && /g/data/va71/crispr_pipeline_annotation/annotation_database/annotation_upgrades/virsorter2/virsort_cmd2.sh " + b + " " + b + "_virsort_all " + str(cores)],shell=True)
	virsorter_parser2.parse(b + "_virsort_all/" + "final-viral-boundary.tsv",b + "_virsort_all/" + "final-viral-score.tsv", b + "_virsort_all/" + "final_table.csv")

	subprocess.run(["source /g/data/va71/my_conda/new_conda/etc/profile.d/conda.sh && conda activate vs2 && /g/data/va71/crispr_pipeline_annotation/annotation_database/annotation_upgrades/virsorter2/virsort_cmd2.sh " + b_phage + " " + b_phage + "_virsort_all " + str(cores)],shell=True)
	virsorter_parser2.parse(b_phage + "_virsort_all/" + "final-viral-boundary.tsv",b_phage + "_virsort_all/" + "final-viral-score.tsv", b_phage + "_virsort_all/" + "final_table.csv")
	exit()

rna_cds_not_run = True
if (rna_cds_not_run):
	subprocess.run(["/g/data/va71/crispr_pipeline_annotation/prodigal/prodigal", "-p", "meta", "-i", b, "-o", b + "_output_full.txt", "-a", b + "_aa_raw.fasta", "-d", b + "_CDS_raw.fnn"]) # output should be redirected to genome specific folder

	
	subprocess.run(["/g/data/va71/crispr_pipeline_annotation/prodigal/prodigal", "-p", "meta", "-i", b_phage, "-o", b_phage + "_output_full.txt", "-a", b_phage + "_aa_raw.fasta", "-d", b_phage + "_CDS_raw.fnn"]) # output should be redirected to genome specific folder
	subprocess.run(["samtools", "faidx", b + "_CDS_raw.fnn"])
	subprocess.run(["samtools", "faidx", b_phage + "_CDS_raw.fnn"])

# END OF RNA CDS PREDICTION
# subprocess.run(["/g/data/va71/crispr_pipeline_annotation/plasme.sh " + b + " " + b + "_plasmid_overlaps.fa " + str(cores)],shell=True)
# subprocess.run(["/g/data/va71/crispr_pipeline_annotation/plasme.sh " + b_phage + " " + b_phage + "_plasmid_overlaps.fa " + str(cores)],shell=True)

my_pickle_phage = open(phage_a + "_pickle_block_dict.pickle", "rb")
phage_block_dict = pickle.load(my_pickle_phage)
my_pickle_phage.close()

if (phage_genomes == 1):
	phage_protein_annotation_frame = []
	phage_summary_frame = phage_summary_frame[1:]
	pool = Pool(cores)
	db_directory_switch = 0
	target_pos_switch = 0



my_pickle = open(a + "_pickle_block_dict.pickle", "rb")
block_dict = pickle.load(my_pickle)
my_pickle.close()
sequence_master_info_file = open(a + "_summary.csv.txt", "r")
summary_frame = list(csv.reader(sequence_master_info_file))
sequence_master_info_file.close()
db_directory_switch = 1

# need to create the table output headers here.

out_table = open (b + "_annotations.csv","a")
spamwriter = csv.writer(out_table)
header_norm = []
header_norm.extend(["genome_id","sequence_id", "protein_start_site", "protein_end_site", "sense","orf_url", "genome_length"])
header_norm.extend(["hhblits_id", "hhblits_description", "hhblits_probability", "hhblits_E-value", "hhblits_score", "hhblits_similarity"])
header_norm.extend(["pfam_blits_id", "pfam_description", "pfam_probability", "pfam_E-value", "pfam_score", "pfam_similarity"])
header_norm.extend(["crisprcasfinder_blits_id", "crisprcasfinder_description", "crisprcasfinder_score", "crisprcasfinder_c_E-value", "crisprcasfinder_i_E-value"])
#header_norm.extend(["padlocplus_blits_id", "padlocplus_description", "padlocplus_score", "padlocplus_c_E-value", "padlocplus_i_E-value"])
header_norm.extend(["type_III_signal_id","type_III_signal_description","type_III_signal_score","type_III_signal_c_E_value","type_III_i_E_value"])
print("Start of HMM prediction:")
spamwriter.writerow(header_norm)
out_table.close()

out_table = open (b_phage + "_annotations.csv","a")
spamwriter = csv.writer(out_table)
header_norm = []
header_norm.extend(["phage_id","sequence_id", "protein_start_site", "protein_end_site", "sense","orf_url", "genome_length"])
header_norm.extend(["hhblits_id", "hhblits_description", "hhblits_probability", "hhblits_E-value", "hhblits_score", "hhblits_similarity"])
header_norm.extend(["pfam_blits_id", "pfam_description", "pfam_probability", "pfam_E-value", "pfam_score", "pfam_similarity"])
header_norm.extend(["crisprcasfinder_blits_id", "crisprcasfinder_description", "crisprcasfinder_score", "crisprcasfinder_c_E-value", "crisprcasfinder_i_E-value"])
#header_norm.extend(["padlocplus_blits_id", "padlocplus_description", "padlocplus_score", "padlocplus_c_E-value", "padlocplus_i_E-value"])
header_norm.extend(["type_III_signal_id","type_III_signal_description","type_III_signal_score","type_III_signal_c_E_value","type_III_i_E_value"])
print("Start of HMM prediction:")
spamwriter.writerow(header_norm)
out_table.close()
out_table = open(b + "_folded.csv", "a")
spamwriter = csv.writer(out_table)
spamwriter.writerow(["Genome_id","Protein_id", "MEA", "Probability"])
out_table.close()
out_table = open(b_phage + "_folded.csv", "a")
spamwriter = csv.writer(out_table)
spamwriter.writerow(["Genome_id","Protein_id", "MEA", "Probability"])
out_table.close()
if (phage_only == 0):
	pool = Pool(cores)
	for frame in summary_frame[1:]:
		pool.apply_async(annotation_parallelisation, (frame, block_dict, db_directory_switch, b, db_directory_path, phmmer_switch,hhsearch_switch,hhblits_switch,dali_switch,target_pos_switch,pfam_switch, criscasfinder_switch, merge_protein,eggnog_switch,vmatch_switch,virsorter_switch,rnafold_switch,padlocplus_switch_novel, False))
	#	pool.apply_async(annotation_parallelisation, (frame, block_dict, db_directory_switch, b, db_directory_path, phmmer_switch,hhsearch_switch,hhblits_switch, dali_switch, target_pos_switch, pfam_switch, criscasfinder_switch, merge_protein, False))
	pool.close()
	pool.join()


if (host_only == 0):
	pool = Pool(cores)
	for frame in phage_summary_frame[1:]:
		pool.apply_async(annotation_parallelisation, (frame, phage_block_dict, db_directory_switch, b_phage, phage_db_directory_path, phmmer_switch,hhsearch_switch,hhblits_switch,dali_switch,target_pos_switch,pfam_switch, criscasfinder_switch, merge_protein,eggnog_switch,vmatch_switch,virsorter_switch,rnafold_switch,padlocplus_switch_novel, True))

	pool.close()
	pool.join()


