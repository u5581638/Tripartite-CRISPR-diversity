
from Bio import SeqIO
import csv
import sys
import os
import subprocess
import phmmer_parser3
import virsorter_parser2
import hhsuite_parser
import dali_parser
import shutil
import timeit
from time import process_time
from time import time
import fnmatch
from multiprocessing import Pool
import pickle
import vmatch_plasmid_detector
import viennafold_parser
import re 
# from multiprocessing import Lock # need something more powerful (Gadi global)
import fcntl
import numpy
import sqlite3
import pandas

def annotation_parallelisation (frame, block_dict, db_directory_switch, b, db_directory_path, is_phage, phmmer_switch=0,hhsearch_switch=0,hhblits_switch=0, dali_switch=0, target_pos_switch=0, pfam_switch=0, criscasfinder_switch=0, merge_protein=0, vmatch_switch=0, eggnog_switch=0, crispr_switch=1,virsorter_switch=0,rnafold_switch=0,alphafold_switch=0, padlocplus_switch_novel=0):
	con = sqlite3.connect("/g/data/va71/crispr_pipeline_annotation/" + "crispr-phage_interactome.sql")
	cur = con.cursor()
	rnafold_switch = 0
	vmatch_switch = 0
	eggnog_switch = 0 
	virsorter_switch = 0 
	hhblits_switch = 1
	criscasfinder_switch = 1
	type_III_signal_switch = 1
	
	if (db_directory_switch != 0):
		print("Hih")
		frame_url_index = 4 # what does this line do? A: It determines which column in the table the genome_id is located!!
		if (is_phage):
			frame_url_index = 3

		folder_name = b + "_genome_annotation_collation_" + frame[frame_url_index]
		folder_name = folder_name.split(".")
		new_folder_name = folder_name[0] + "_" + folder_name[-1]
		folder_name = new_folder_name + "/"

		if (os.path.isdir(folder_name)):
			pass 
		else:	
			os.mkdir(folder_name)
		#	print(folder_name)
		#	exit()

		# this is the line that fails!!!!
		print("Progress 1")
	#	print(frame)
		if frame[frame_url_index] in block_dict:
			print("Match 1!")
			genome = block_dict[frame[frame_url_index]]
			genome_len = len(genome)
			protein_file_name = folder_name  + "/" + frame[frame_url_index] + ".fasta" 
			SeqIO.write(genome, protein_file_name, "fasta")

		else:
			print ("error, genome key not recognised")
			print(frame[frame_url_index])
			print("Whole frame:")
			print(frame)
	else:
		frame_url_index = 3
		folder_name = b + "_genome_annotation_collation_" + frame[frame_url_index]
		folder_name = folder_name.split(".")
		new_folder_name = folder_name[0] + "_" + folder_name[-1]

		folder_name = new_folder_name + "/"
		if (os.path.isdir(folder_name)):
			pass 
		else:	
		#	print("folder_dir:")
		#	print(folder_name)
			os.mkdir(folder_name)
		if frame[frame_url_index] in block_dict:
			genome = block_dict[frame[frame_url_index]]
			genome_len = len(genome)
			protein_file_name = folder_name  + "/" + frame[frame_url_index] + ".fasta" # calling this variable "protein file name is confusing!!!"
			SeqIO.write(genome, protein_file_name, "fasta")
			print("protein_created!!")
		else:
			print ("error, genome key not recognised")
			print(frame)


	protein_annotation_frame = []

	# START OF PHAGE AND GENERAL GENOME PREDICTIONS COUPLED TO ANNOTATION!!

	# VMATCH -> circular DNA detection. Priority!! Might be the hardest capability to integrate
	if (vmatch_switch == 1):
		print(folder_name)
		print(protein_file_name)
		my_dir = os.getcwd()
		os.chdir(folder_name)
		subprocess.run(["/g/data/va71/vmatch-2.3.1-Linux_x86_64-64bit/mkvtree", "-db", protein_file_name,"-pl","-allout", "-dna", "-v"])
		subprocess.run(["/g/data/va71/vmatch-2.3.1-Linux_x86_64-64bit/vmatch" + " -p" + " -l" + " 10" + " -d " + protein_file_name + " > " + protein_file_name + "_vmatch_drs.txt"],shell=True)
		subprocess.run(["rm *.al1 *.bck *.bwt *.des *.lcp *.llv *.ois *.prj *.sds *.skp *.stil *.suf *.tis"],shell=True)
		os.chdir(my_dir)
		my_arg1 = protein_file_name + "_vmatch_drs.txt"
		my_arg2 = protein_file_name
		my_arg3 = protein_file_name + "_direct_inverted_pallindromic_repeats.csv"
		vmatch_plasmid_detector.plasmid_detect(my_arg1, my_arg2, my_arg3) # This needs to generate a reduced list of circular and linear genomes, as well as a list of direct and inverted repeats in csv format. This should be saved in the same directory (genomes or spacer_distribution_analysis)
		vmatch_url = open(protein_file_name + "_direct_inverted_pallindromic_repeats.csv"+ "_drs.csv","r")
		vmatch_table = pandas.read_csv(vmatch_url)

		fcntl.flock(vmatch_url.fileno(), fcntl.LOCK_EX)
		print("Check my phages!!!!!")
		print(is_phage)
		'''
		if (is_phage):
			vmatch_table.rename(columns={"Genome_id":"Phage_id"},inplace=True)
			# need to add a command 

			cur.execute("CREATE TABLE IF NOT EXISTS VMATCH_REPEATS_PHAGE (Phage_id STRING, length STRING, sequence_number STRING, relative_position STRING, type STRING, length2 STRING, sequence_number2 STRING, relative_position2 STRING, distance_value STRING, `E-value` STRING, score STRING, `perc-identity` STRING)")
			genome_exists = cur.execute("SELECT PHAGE_ID FROM VMATCH_REPEATS_PHAGE WHERE PHAGE_ID=?",(vmatch_table.iloc[0][0],))
			genome_row = genome_exists.fetchone()
			if (genome_row is None):
				vmatch_table.to_sql('VMATCH_REPEATS_PHAGE',con, if_exists='append', index=False)
		else:	
			cur.execute("CREATE TABLE IF NOT EXISTS VMATCH_REPEATS (Genome_id STRING, length STRING, sequence_number STRING, relative_position STRING, type STRING, length2 STRING, sequence_number2 STRING, relative_position2 STRING, distance_value STRING, `E-value` STRING, score STRING, `perc-identity` STRING)")
			genome_exists = cur.execute("SELECT GENOME_ID FROM VMATCH_REPEATS WHERE GENOME_ID=?",(vmatch_table.iloc[0][0],))
			genome_row = genome_exists.fetchone()
			if (genome_row is None):
				vmatch_table.to_sql('VMATCH_REPEATS',con, if_exists='append', index=False)
		'''
		con.commit()
		 # use protein file name instead of * if * doesn't work!!
		vmatch_url.close()
	# Add vmatch results to sql. Ensure output csv file includes genome_id!!
	print("Got to here!!!")
	# EGGNOG2: THIS MIGHT BE THE EASIEST CAPABILITY TO INTEGRATE (IF OUTPUT IS EXCEL!!)
	# NOTE: EGGNOG not used.
	if (eggnog_switch == 1):
		# Eggnog needs to be pre-installed. It seems that only a very particular combination of flags works!!
		# Eggnog may unfortunately be two slow to be useful!!
		# Good way to cross-validate predictions
		subprocess.run(["emapper.py", "-m", "diamond", "--itype", "CDS", "-i", protein_file_name, "-o", protein_file_name + "_eggnog", "--excel", "--sensmode", "fast", "--data_dir", "/g/data/va71/eggnog/eggnog/data/"])
		eggnog_table_url = open(protein_file_name + "_eggnog.emapper.annotations.xlsx","rb")
		eggnog_table = pandas.read_excel(eggnog_table_url, sheet_name='Sheet1')
		eggnog_table = eggnog_table.iloc[2:-3]
		eggnog_table = eggnog_table.set_axis(["query", "seed_ortholog", "evalue","score", "eggNOG_OGs", "max_annot_lvl", "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy","BiGG_Reaction", "PFAMs"], axis=1)
		eggnog_table.rename(columns={"query":"Genome_id"},inplace=True)
		print("Check 3 dataframe")

		fcntl.flock(eggnog_table_url.fileno(), fcntl.LOCK_EX)
		eggnog_table.to_csv(protein_file_name + "_eggnog.emapper.annotations.csv")
		con.commit()
		eggnog_table_url.close()
		# upload excel output to sql via pandas!!

	# virsorter prediction.
	if (virsorter_switch == 1):
		# unsure if this is sufficent?????
		subprocess.run(["source /g/data/va71/my_conda/new_conda/etc/profile.d/conda.sh && conda activate vs2 && /g/data/va71/crispr_pipeline_annotation/annotation_database/annotation_upgrades/virsorter2/virsort_cmd.sh " + protein_file_name + " " + protein_file_name + "_virsort", 1],shell=True)
		virsorter_parser2.parse(protein_file_name + "_virsort/" + "final-viral-boundary.tsv",protein_file_name + "_virsort/" + "final-viral-score.tsv", protein_file_name + "_virsort/" + "final_table.csv")
		virsorter_url = open(protein_file_name + "_virsort/" + "final_table.csv","r")
		
		fcntl.flock(virsorter_url.fileno(),fcntl.LOCK_EX)
		'''
		virsorter_table = pandas.read_csv(virsorter_url)
		if (is_phage):
			virsorter_table.rename(columns={"seqname":"Phage_id"},inplace=True)
			virsorter_table.drop(columns=virsorter_table.columns[0], axis=1, inplace=True)
			cur.execute("CREATE TABLE IF NOT EXISTS VIRSORTER_PREDICTIONS_PHAGE (Phage_id STRING PRIMARY KEY, trim_orf_index_start STRING, trim_orf_index_end STRING,trim_bp_start STRING,trim_bp_end STRING,trim_pr STRING, trim_pr_max STRING, prox_orf_index_start STRING, prox_orf_index_end STRING,prox_bp_start STRING,prox_bp_end STRING,prox_pr STRING, prox_pr_max STRING,partial STRING,full_orf_index_start STRING,full_orf_index_end STRING,full_bp_start STRING,full_bp_end STRING,pr_full STRING,arc STRING,bac STRING,euk STRING,vir STRING, mix STRING, unaligned STRING, hallmark_cnt STRING, \"group\" STRING, shape STRING, seqname_new STRING, final_max_score STRING, final_max_score_group STRING,dsDNAphage STRING,ssDNA STRING, FOREIGN KEY (Phage_id) REFERENCES PHAGE_GENOMES(Phage_id))")
			if (len(virsorter_table["Phage_id"]) != 0):
				genome_exists = cur.execute("SELECT PHAGE_ID FROM VIRSORTER_PREDICTIONS_PHAGE WHERE PHAGE_ID=?",(virsorter_table.iloc[0][0],))
				genome_row = genome_exists.fetchone()
				if (genome_row is None):
					virsorter_table.to_sql("VIRSORTER_PREDICTIONS_PHAGE", con, if_exists='append',index=False)
		else:
			virsorter_table.rename(columns={"seqname":"Genome_id"},inplace=True)
			virsorter_table.drop(columns=virsorter_table.columns[0], axis=1, inplace=True)	
			cur.execute("CREATE TABLE IF NOT EXISTS VIRSORTER_PREDICTIONS (Genome_id STRING PRIMARY KEY, trim_orf_index_start STRING, trim_orf_index_end STRING,trim_bp_start STRING,trim_bp_end STRING,trim_pr STRING, trim_pr_max STRING, prox_orf_index_start STRING, prox_orf_index_end STRING,prox_bp_start STRING,prox_bp_end STRING,prox_pr STRING, prox_pr_max STRING,partial STRING,full_orf_index_start STRING,full_orf_index_end STRING,full_bp_start STRING,full_bp_end STRING,pr_full STRING,arc STRING,bac STRING,euk STRING,vir STRING, mix STRING, unaligned STRING, hallmark_cnt STRING, \"group\" STRING, shape STRING, seqname_new STRING, final_max_score STRING, final_max_score_group STRING,dsDNAphage STRING,ssDNA STRING, FOREIGN KEY (Genome_id) REFERENCES GENOMES(Genome_id))")
			if (len(virsorter_table["Genome_id"]) != 0):
				genome_exists = cur.execute("SELECT GENOME_ID FROM VIRSORTER_PREDICTIONS_PHAGE WHERE GENOME_ID=?",(virsorter_table.iloc[0][0],))
				genome_row = genome_exists.fetchone()
				if (genome_row is None):
					virsorter_table.to_sql("VIRSORTER_PREDICTIONS", con, if_exists='append',index=False) # If appending, it may accept the initialised genome_id file name
		
		con.commit()
		'''
		virsorter_url.close()
	#	subprocess.run(["rm -r " + protein_file_name + "_virsort/"],shell=True)

	# protein_file_name_1 is the name of the genome. 
	# To get protein name ->
	# open file.
	# use re to return matching lines
	if (os.path.isfile(protein_file_name + "_aa_raw.fasta")):
		os.remove(protein_file_name + "_aa_raw.fasta")
	if (os.path.isfile(protein_file_name + "_CDS_raw.fnn")):
		os.remove(protein_file_name + "_CDS_raw.fnn")

	my_fai = open(b + "_aa_raw.fasta" + ".fai", "r")
	fai_read = my_fai.readlines() # read with \n seperated
	headlines = [] 
	for line in fai_read:
		my_id = genome.id
		my_id = genome.id.split("|")
		if (len(my_id) == 2):
			my_id = my_id[0] + '\|' + my_id[1]
		else:
			my_id = my_id[0]	  

		headline = re.search(my_id + '.*', line)

		if (headline is not None):
			print("Yaya!!")
			my_headline = headline.group(0)
			my_headline = my_headline.split("\t") [0]

			subprocess.run(["/g/data/va71/crispr_pipeline_annotation/seqkit faidx " + b + "_aa_raw.fasta " + "\'" + my_headline + "\'" + " -f " + ">> " + protein_file_name + "_aa_raw.fasta"], shell=True) # will only work if the samtools code is italised

	# index prodigal predicted CDS sequences
	if (rnafold_switch == 1): 
		my_fai_cds = open(b + "_CDS_raw.fnn" + ".fai", "r")
		fai_read = my_fai_cds.readlines() # read with \n seperated
		headlines = [] 
		for line in fai_read:
			my_id = genome.id
			my_id = genome.id.split("|")
			if (len(my_id) == 2):
				my_id = my_id[0] + '\|' + my_id[1]
			else:
				my_id = my_id[0]	  

			headline = re.search(my_id + '.*', line)
			if (headline is not None):
				my_headline = headline.group(0)
				my_headline = my_headline.split("\t") [0]
				subprocess.run(["/g/data/va71/crispr_pipeline_annotation/seqkit faidx " + b + "_CDS_raw.fnn " + "\'" + my_headline + "\'" + " -f " + ">> " + protein_file_name + "_CDS_raw.fnn"], shell=True) # will only work if the samtools code is italised
		
	print("finished the outer loop part of the protein annotation!!!")


	my_fai.close()

	# main code here for gene annotation tabulation
	print("run_prodigal!!")
	# before the file is opened in append mode, the program should check whether the file exists.
	# This should be true for all file initiation in append mode!!
	if (os.path.isfile(b + "_genome_annotation_collation/" + frame[frame_url_index])):
		os.remove(b + "_genome_annotation_collation/" + frame[frame_url_index])
	if (os.path.isfile(protein_file_name + "_aa_raw.fasta")):
		predicted_proteins = list(SeqIO.parse(protein_file_name + "_aa_raw.fasta", "fasta"))
	else:
		print("Error: No protein file found!!")	
	if (rnafold_switch == 1):
		cds_proteins = list(SeqIO.parse(protein_file_name + "_CDS_raw.fnn", "fasta"))

	# this table will need a header!!!
	# check this section of code!!
	if (os.path.isdir(folder_name + "hmmer_sequence_files/")):
		pass 
	else:	
		hmm_file_dir = os.mkdir(folder_name + "hmmer_sequence_files/")
	if (alphafold_switch == 1):
		if (os.path.isdir(db_directory_path + "protein_predictions/")):
			pass 
		else:
			alpha_file_dir = os.mkdir(db_directory_path + "protein_predictions/")
	if (rnafold_switch == 1):
		if (os.path.isdir(db_directory_path + "RNAfold_predictions/")):
			print("I passed!")
			print(db_directory_path)
			pass 
		else:
			print("made_dir_of_some_description")
			print(db_directory_path)
			rnafold_dir = os.mkdir(db_directory_path + "RNAfold_predictions/")		

	if (is_phage):
		protein_annotation_frame.append(["phage_id","sequence_id", "protein_start_site", "protein_end_site", "sense","orf_url", "genome_length"])
	else:	
		protein_annotation_frame.append(["genome_id","sequence_id", "protein_start_site", "protein_end_site", "sense","orf_url", "genome_length"])
	
	if (phmmer_switch == 1):
		protein_annotation_frame[0].extend(["phmmer_prediction_protein_ID", "phmmer_prediction", "phmmer_score", "c_e-value", "i_e_value"])
	if (hhsearch_switch == 1):
		protein_annotation_frame[0].extend(["hhsearch_id", "hhsearch_description", "probability", "E-value", "score", "similarity"])
	if (hhblits_switch == 1):
		protein_annotation_frame[0].extend(["hhblits_id", "hhblits_description", "hhblits_probability", "hhblits_E-value", "hhblits_score", "hhblits_similarity"])
	if (dali_switch == 1):
		protein_annotation_frame[0].extend(["AF2+DALI prediction"])
	if (target_pos_switch == 1):
		protein_annotation_frame[0].extend(["target_protein_start", "target_protein_end", "target_protein_sense", "target_protein_order"])	
	if (pfam_switch == 1):
		protein_annotation_frame[0].extend(["pfam_blits_id", "pfam_description", "pfam_probability", "pfam_E-value", "pfam_score", "pfam_similarity"])	
	if (criscasfinder_switch == 1):
		protein_annotation_frame[0].extend(["crisprcasfinder_blits_id", "crisprcasfinder_description", "crisprcasfinder_score", "crisprcasfinder_c_E-value", "crisprcasfinder_i_E-value"])
	if (padlocplus_switch_novel == 1):
		protein_annotation_frame[0].extend(["padlocplus_blits_id", "padlocplus_description", "padlocplus_score", "padlocplus_c_E-value", "padlocplus_i_E-value"])	
	if (type_III_signal_switch == 1):
		protein_annotation_frame[0].extend(["type_III_signal_id","type_III_signal_description","type_III_signal_score","type_III_signal_c_E_value","type_III_i_E_value"])
	print("Start of HMM prediction:")
	print(len(predicted_proteins))
	print(frame)


	a = 0
	while (a < len(predicted_proteins)):
		sequence_description = predicted_proteins[a].description.split(" # ")
		sequence_identifer = sequence_description[4]
		sequence_identifer_components = sequence_identifer.split(";")
		sequence_identifer = sequence_identifer_components[0]
		protein_start_site = sequence_description[1]
		protein_end_site = sequence_description[2]
		protein_sense = sequence_description[3]
		fourth_frame = frame[frame_url_index]
		sequence_identifer = fourth_frame + "_" + sequence_identifer.split("ID=") [1]
		first_frame = frame[1]
		
		protein_annotation_frame.append([fourth_frame,sequence_identifer,protein_start_site,protein_end_site,protein_sense, first_frame,str(genome_len)])
		

		# section to run HHpred and HMMER within the same loop!!

		sequence_id = predicted_proteins[a].id 
		sequence_id = sequence_id.split("|")
		sequence_id = sequence_id[-1] # splitting this might cause problems with sql!!
		SeqIO.write(predicted_proteins[a], folder_name + "hmmer_sequence_files/" + sequence_id, "fasta")
		# Write CDS genome sequences to file
		if (rnafold_switch == 1):
			SeqIO.write(cds_proteins[a], db_directory_path + "RNAfold_predictions/" + sequence_id, "fasta")
		# RNA secondary structure prediction of ORFs using viennafold
		if (rnafold_switch == 1):
			subprocess.run(["/g/data/va71/crispr_pipeline_annotation/viennafold/ViennaRNA-2.5.1/src/bin/RNAfold " + "-p " + "--MEA " + "-i " + folder_name + "hmmer_sequence_files/" + sequence_id + " > " + db_directory_path + "RNAfold_predictions/" + sequence_id + "folded.txt"], shell=True)
			viennafold_parser.parse(db_directory_path + "RNAfold_predictions/" + sequence_id + "folded.txt", annotations_prefix + "_folded.csv")
			os.remove(db_directory_path + "RNAfold_predictions/" + sequence_id)
			os.remove(db_directory_path + "RNAfold_predictions/" + sequence_id + "folded.txt")
			subprocess.run(["rm /g/data/va71/crispr_pipeline_annotation/annotation_upgraded_main_workflow_run/running_scripts/annotation_only/*.ps"],shell=True)
			# add statement to add RNAfold results to Sql here!!
		#	my_vienna = open(db_directory_path + "RNAfold_predictions/" + sequence_id + "folded.csv","r")
		#	vienna_sql = pandas.read_csv(my_vienna)
		#	fcntl.flock(my_vienna.fileno(), fcntl.LOCK_EX) # This has been verified to work)
			'''
			if (is_phage):
				cur.execute('CREATE TABLE IF NOT EXISTS RNAFOLD_PREDICTIONS_PHAGE (Protein_id STRING,Hairpin_number STRING, Minimum_Free_Energy STRING)')
				vienna_sql.to_sql('RNAFOLD_PREDICTIONS_PHAGE', con, if_exists='append',index=False)
			else: 
				cur.execute('CREATE TABLE IF NOT EXISTS RNAFOLD_PREDICTIONS (Protein_id STRING,Hairpin_number STRING, Minimum_Free_Energy STRING)')
				vienna_sql.to_sql('RNAFOLD_PREDICTIONS', con, if_exists='append',index=False)
			'''
		#	my_vienna.close()

		# HMM prediction using phmmer search against uniref90 database.
		if (phmmer_switch == 1):
			t1_start = time()
			subprocess.run(["/g/data/va71/alphafold/install/hmmer-3.3.2/src/phmmer", "-o", folder_name + "hmmer_sequence_files/" + sequence_id + "_phmmer_hits.fa", "--cpu", "1", folder_name + "hmmer_sequence_files/" + sequence_id, "/g/data/va71/alphafold/datasets/alphafold2/uniref90/uniref90.fasta"] )
			t1_end = time()
			print("time to run phmmer is:", t1_start, t1_end, t1_end - t1_start)
			phmmer_frames = phmmer_parser3.phmmer_parse(folder_name + "hmmer_sequence_files/" + sequence_id + "_phmmer_hits.fa") #only works for uniref90 as of 10/11/2021
			phmmer_frames_zero = phmmer_frames[0][0]
			phmmer_frames_one = phmmer_frames [0][1]
			score = phmmer_frames[0][2]
			c_e_value = phmmer_frames[0][3]
			i_e_value = phmmer_frames[0][4]
			protein_annotation_frame[-1].extend([phmmer_frames_zero, phmmer_frames_one, score, c_e_value, i_e_value])
			print("phmmer done!!")
		# HMM prediction using hhsearch against pdb70 database.
		if (hhsearch_switch == 1):
		#	protein_annotation_frame[0].extend(["hhsearch_id", "hhsearch_description"])
			t1_start = time()
			subprocess.run(["/g/data/va71/alphafold/install/hh-suite/build/src/hhsearch", "-i", folder_name + "hmmer_sequence_files/" + sequence_id, "-o",  folder_name + "hmmer_sequence_files/" + sequence_id + "_hhsearch_hits.txt.fa", "-cpu", "1", "-e", "0.001", "-maxseq", "5000", "-d", "/g/data/va71/alphafold/datasets/alphafold2/pdb70/pdb70"])
			t1_end = time()
			print("time to run hhsearch is", t1_start, t1_end, t1_end - t1_start)
			hhsearch_frames = hhsuite_parser.pdb70_parse(folder_name + "hmmer_sequence_files/" + sequence_id + "_hhsearch_hits.txt.fa") # this will only work for pdb70/pdb70. Need to write this script tonight!!
			hhsearch_frames_zero = hhsearch_frames[0][0]
			hhsearch_frames_one = hhsearch_frames[0][1]
		#	print(hhsearch_frames)
			hhsearch_frames_probability = hhsearch_frames[0][2]
			hhsearch_frames_e_value = hhsearch_frames[0][3]
			hhsearch_frames_score = hhsearch_frames[0][4] 
			hhsearch_frames_similarity = hhsearch_frames[0][5] 
			protein_annotation_frame[-1].extend([hhsearch_frames_zero, hhsearch_frames_one, hhsearch_frames_probability, hhsearch_frames_e_value, hhsearch_frames_score])
			print("hhsearch complete!!")
		# HMM prediction using HHBlits against the Pdb70 database
		if (hhblits_switch == 1):
		#	protein_annotation_frame[0].extend(["hhblits_id", "hhblits_description"])
			t1_start = time()
			subprocess.run(["/g/data/va71/alphafold/install/hh-suite/build/src/hhblits", "-i", folder_name + "hmmer_sequence_files/" + sequence_id, "-o",  folder_name + "hmmer_sequence_files/" + sequence_id + "_hhblits_hits.txt.fa", "-cpu", "1", "-e", "0.001", "-maxseq", "15000", "-d", "/g/data/va71/alphafold/datasets/alphafold2/pdb70/pdb70", "-v", "0"])
			t1_end = time()
			print("time to run hhblits is:", t1_start, t1_end, t1_start - t1_end)
			print("hhblits complete!!")	
			hhblits_frames = hhsuite_parser.pdb70_parse(folder_name + "hmmer_sequence_files/" + sequence_id + "_hhblits_hits.txt.fa") # this will only work for pdb70. Need to write this script tonight!! Will use biopython!!
			if (hhblits_frames == []):

				hhblits_frames_zero = "NA"
				hhblits_frames_one = "NA"
				hhblits_frames_probability = "NA"
				hhblits_frames_e_value = "NA"
				hhblits_frames_score = "NA"
				hhblits_frames_similarity = "NA"
				protein_annotation_frame[-1].extend([hhblits_frames_zero, hhblits_frames_one, hhblits_frames_probability, hhblits_frames_e_value, hhblits_frames_score, hhblits_frames_similarity])
			else:			
				hhblits_frames_zero = hhblits_frames[0][0]
				hhblits_frames_one = hhblits_frames[0][1]
				hhblits_frames_probability = hhblits_frames[0][2]
				hhblits_frames_e_value = hhblits_frames[0][3]
				hhblits_frames_score = hhblits_frames[0][4]
				hhblits_frames_similarity = hhblits_frames[0][5]
				protein_annotation_frame[-1].extend([hhblits_frames_zero, hhblits_frames_one, hhblits_frames_probability, hhblits_frames_e_value, hhblits_frames_score, hhblits_frames_similarity])
		# run alphafold followed by a dali structure-based homology search. Note: This feature did not end up being implemented.
		if (dali_switch == 1):
		#	protein_annotation_frame[0].extend(["AF2+DALI prediction"])
			t1_start = time()
			subprocess.run(["./alphafold_gpu_enabled_cpu_dali.sh",folder_name + "hmmer_sequence_files/" + sequence_id ]) # run_alphafold
			# now run dali
			subprocess.run(["DaliLite.v5/bin/import.pl", "--pdbfile", sequence_id + "/" + "ranked_0.pdb"])
			subprocess.run(["DaliLite.v5/bin/dali.pl", "--dat1", "./DAT/","--dat2", "./DAT/"]) ## need to download a database to run this line!!!
			dali_prediction = dali_parser.parse("#output_table_from_dali#!!!")
			protein_annotation_frame[-1].extend(dali_prediction)
		'''
		if (target_pos_switch == 1):
			protein_annotation_frame[-1].extend([genome_target_start, genome_target_end, genome_target_sense, genome_target_order])	
		'''

		# HMM prediction using HHBlits against the Pfam database
		if (pfam_switch == 1):
			t1_start = time()
			subprocess.run(["/g/data/va71/alphafold/install/hh-suite/build/src/hhblits", "-i", folder_name + "hmmer_sequence_files/" + sequence_id, "-o",  folder_name + "hmmer_sequence_files/" + sequence_id + "_hhblits_hits.txt.fa", "-cpu", "1", "-e", "0.001", "-maxseq", "15000", "-d", "/g/data/va71/alphafold/datasets/alphafold2/pfam2/pfam", "-v", "0"])
			t1_end = time()
			print("time to run hhblits pfam is:", t1_start, t1_end, t1_start - t1_end)
			print("hhblits complete!!")	
		# parsers go here.
		# need to include HHsearch + HHblits parsers to complete protein annotation
			hhblits_frames = hhsuite_parser.pdb70_parse_pfam(folder_name + "hmmer_sequence_files/" + sequence_id + "_hhblits_hits.txt.fa") # this will only work for pdb70. Need to write this script tonight!! Will use biopython!!

			if (hhblits_frames == []):

				hhblits_frames_zero = "NA"
				hhblits_frames_one = "NA"
				hhblits_frames_probability = "NA"
				hhblits_frames_e_value = "NA"
				hhblits_frames_score = "NA"
				hhblits_frames_similarity = "NA"
				protein_annotation_frame[-1].extend([hhblits_frames_zero, hhblits_frames_one, hhblits_frames_probability, hhblits_frames_e_value, hhblits_frames_score, hhblits_frames_similarity])
			else:			
				hhblits_frames_zero = hhblits_frames[0][0]
				hhblits_frames_one = hhblits_frames[0][1]
				hhblits_frames_probability = hhblits_frames[0][2]
				hhblits_frames_e_value = hhblits_frames[0][3]
				hhblits_frames_score = hhblits_frames[0][4]
				hhblits_frames_similarity = hhblits_frames[0][5]
				protein_annotation_frame[-1].extend([hhblits_frames_zero, hhblits_frames_one, hhblits_frames_probability, hhblits_frames_e_value, hhblits_frames_score, hhblits_frames_similarity])
		
		# add code here to run code tabulating padloc annotations using hmmscan. change url directories to padloc database. Add the ability to add new HMMs. This must cause all running threads to temporarily pause until the HMM addition has finished writing.
		
		# In practice this should only work for CRISPR-associated proteins - NOT phage proteins unless a reliable database of 
		# add new HMM automatically for unknown ORFs. Need to iteratively build HMM from one sequence but augment with others using 1 sequence HMM or BLAST.. Need to declare a new HMM file.
		# Note, this functionality did not end up being used.
		if (padlocplus_switch_novel == 1): 
			
			t1_start = time()
			subprocess.run(["/g/data/va71/alphafold/install/hmmer-3.3.2/src/hmmscan", "-o", folder_name + "hmmer_sequence_files/" + sequence_id + "_hmmscan_hits.fa", "--cpu", "1", "/g/data/va71/alphafold/datasets/alphafold2/crisprcasfinder_hmm_profiles/non_dup_merged_all_models_master_TIGR.HMM", folder_name + "hmmer_sequence_files/" + sequence_id] )
			t1_end = time()
			print("time to run hmmscan is:", t1_start, t1_end, t1_end - t1_start)
			# phmmer_parser3
			phmmer_frames = hmm_scan_parser.hmmscan_parse(folder_name + "hmmer_sequence_files/" + sequence_id + "_hmmscan_hits.fa") #only works for uniref90 as of 10/11/2021
		#	may need 2nd hmmscan cmd for novel assigned proteins.	
			print("Working fine!!")
			print(phmmer_frames)
			if (phmmer_frames == []):
				print("Last stand!!")
				phmmer_frames_zero = "NA"
				phmmer_frames_one = "NA"
				score = "NA"
				c_e_value = "NA"
				i_e_value = "NA"
				
				if (phmmer_frames_zero == "NA"): 

				#	current_index =  retrieve from SQL. This step is no longer nessessary I think???
					# steps:
					# 1. Isolate 1 sequence HMM.
					novel_sequence_id = folder_name + "hmmer_sequence_files/" + sequence_id

					# 2. Assign a tmp name. - may not be nessessary!

					# 3. BLAST/Hmmscan

					subprocess.run(["blastp", "-db",b + "_aa_raw.fasta", "-query", novel_sequence_id, "-outfmt", "10", "-out", novel_sequence_id + "_novel_hits.csv", "-max_target_seqs", "200000000", "-max_hsps", "1", "-evalue", "0.0000001" ])
					# 4. Assign a permanent auto generated name if a certain threshold of non-redundant (convert NCBI to JGI numbers using GOLD) hits + 2 sample numbers is achieved (this criteron may be customisable)!
					# This function should retrieve the proteins from the main host/phage protein file
					
					
					with open(novel_sequence_id + "_novel_hits.csv", "r") as csvfile:
						hit_table = list(csv.reader(csvfile)) # load protein hits
					# run seqkit to extract the protein sequences from the first row
					novel_sequence_id = SeqIO.read(novel_sequence_id, "fasta")
					for hit in hit_table:
						print(hit[1])
						# may want to deduplicate. This step is very compute expensive!!!!!!!!
					#	subprocess.run(["/g/data/va71/crispr_pipeline_annotation/seqkit faidx " + b + "_aa_raw.fasta" + " " + "\'" + hit[1] + "\'" + " -f " + ">> " + db_directory_path + "\'" + novel_sequence_id.id + "\'" ],shell=True) # "\'" + novel_sequence_id.id + "\'" + "_novel_cas.fasta"
					
					proteins_out_url = db_directory_path + novel_sequence_id.id
					protein_set = list(SeqIO.parse(proteins_out_url, "fasta"))
					# This issue is in this line!!
					if (len(protein_set) > 1):	
						protein_ga_ids = set()
						for protein in protein_set:
							prot_id = protein.id.split("|")
							if (len(prot_id) > 1): # This may mess up this part of the tool up!!!! COME BACK TO THIS!!
								prot_id = prot_id [1]
								if (prot_id[:2] == "Ga"):
									my_id = prot_id.split("_") [0]
									if my_id not in protein_ga_ids:
										protein_ga_ids.add(my_id)
		
						if (len(protein_ga_ids) > 1):
									

						
					# 5. MSA generation
							subprocess.run(["/g/data/va71/crispr_pipeline_annotation/clustalo -i " + db_directory_path + "\'" + novel_sequence_id.id + "\'" + " > " + db_directory_path + novel_sequence_id.id.split("|")[1]  + "_aligned.fasta"], shell=True) #  
					# 6. Build hmm
							permanent_index_url = open("/g/data/va71/crispr_pipeline_annotation/permanent_index_file.txt", "r")
							index = str(permanent_index_url.read().strip())
							permanent_index_url.close()
							new_index = str(int(index) + 1)
							write_index = str(new_index)
							# need to lock other threads from interfering in this step!!
						#	lock = Lock()
						#	lock.acquire()
							permanent_index_url = open("/g/data/va71/crispr_pipeline_annotation/permanent_index_file.txt","w") # may want to store this on SQL!!!!!
							fcntl.flock(permanent_index_url.fileno(), fcntl.LOCK_EX)
							permanent_index_url.write(write_index) # this can be made doubly safe by restricting the amount written to < 4096 bytes (using a loop)!!!
							permanent_index_url.close() # This should release the lock. Hopefully this lock affects every Gadi process, not just individual nodes. Otherwise have to use MPI!
						#	lock.release()
							in_arg = db_directory_path + novel_sequence_id.id.split("|")[1]  + "_aligned.fasta"
							out_arg = db_directory_path + novel_sequence_id.id.split("|")[1] +  "new_hmm.hmm "
							subprocess.run(["/g/data/va71/alphafold/install/hmmer-3.3.2/src/hmmbuild " + out_arg + in_arg], shell=True)



					# use and update current index at this step (to assign a unique name to the hmm i.e. GAETAN1 ... GAETAN2 )

					# 7. lock other processes
					# 8. Concate to primary padlocdb.hmm file
						#	lock2 = Lock()
						#	lock2.acquire()
							novel_hmm_url = open ("/g/data/va71/crispr_pipeline_annotation/novel_crispr_associated_profiles.hmm", "a")
							fcntl.flock(novel_hmm_url.fileno(), fcntl.LOCK_EX) # This has been verified to work!!
							new_hmm = open(db_directory_path + novel_sequence_id.id.split("|")[1] + "new_hmm.hmm", "r")
							hmm_str = new_hmm.read() + "\n" # may not need the end of line character!!
							novel_hmm_url.write(hmm_str) # this can be made doubly safe by restricting the amount written to < 4096 bytes (using a loop)!!!
							new_hmm.close()
						#	subprocess.run(["cat " + "new_hmm.hmm" + " >> " + "/g/data/va71/crispr_pipeline_annotation/novel_crispr_associated_profiles.hmm"], shell=True) # need to change padlocdb.hmm to actual name of hmm file!!! May consist of only novel hmm files!!
							subprocess.run(["rm -r " + in_arg + " " + out_arg + " " + db_directory_path + "\'" + novel_sequence_id.id + "\'"], shell=True)

					# 9.

					# Allow other proccesses to resume.
						#	lock2.release()
							novel_hmm_url.close()
			else:
				phmmer_frames_zero = phmmer_frames[0][0]
				phmmer_frames_one = phmmer_frames [0][1]
				score = phmmer_frames[0][2]
				c_e_value = phmmer_frames[0][3]
				i_e_value = phmmer_frames[0][4]
			protein_annotation_frame[-1].extend([phmmer_frames_zero, phmmer_frames_one, score, c_e_value, i_e_value])
			# end of annotation


		# gene annotations using HMMscan against the DEFLOC database.
		if (criscasfinder_switch == 1 or (padlocplus_switch_novel == 1 and is_phage)):
			t1_start = time()
			subprocess.run(["/g/data/va71/alphafold/install/hmmer-3.3.2/src/hmmscan", "-o", folder_name + "hmmer_sequence_files/" + sequence_id + "_hmmscan_hits.fa", "--cpu", "1", "/g/data/va71/alphafold/datasets/alphafold2/crisprcasfinder_hmm_profiles/non_dup_merged_all_models_master_TIGR.HMM", folder_name + "hmmer_sequence_files/" + sequence_id] )
			t1_end = time()
			print("time to run hmmscan is:", t1_start, t1_end, t1_end - t1_start)
			# phmmer_parser3
			phmmer_frames = hmm_scan_parser.hmmscan_parse(folder_name + "hmmer_sequence_files/" + sequence_id + "_hmmscan_hits.fa") #only works for uniref90 as of 10/11/2021
			if (phmmer_frames == []):
				phmmer_frames_zero = "NA"
				phmmer_frames_one = "NA"
				score = "NA"
				c_e_value = "NA"
				i_e_value = "NA"
				protein_annotation_frame[-1].extend([phmmer_frames_zero, phmmer_frames_one, score, c_e_value, i_e_value])
			else:
				phmmer_frames_zero = phmmer_frames[0][0]
				phmmer_frames_one = phmmer_frames [0][1]
				score = phmmer_frames[0][2]
				c_e_value = phmmer_frames[0][3]
				i_e_value = phmmer_frames[0][4]
				protein_annotation_frame[-1].extend([phmmer_frames_zero, phmmer_frames_one, score, c_e_value, i_e_value])
		# Gene annotations using HMMscan against HMM profiles containing the cyclase domain from type III CRISPR-Cas systems. This was for personal curiosity rather than the whole thesis.
		if (type_III_signal_switch == 1):
			t1_start = time()
			subprocess.run(["/g/data/va71/alphafold/install/hmmer-3.3.2/src/hmmscan", "-o", folder_name + "hmmer_sequence_files/" + sequence_id + "_hmmscan_hits.fa", "--cpu", "1", "/g/data/va71/crispr_pipeline_annotation/type_III_signalling_profile_construction/ggdef_coA_signal_domain_all.hmm", folder_name + "hmmer_sequence_files/" + sequence_id] )
			t1_end = time()
			print("time to run typeIII hmmscan is:", t1_start, t1_end, t1_end - t1_start)
			# phmmer_parser3
			phmmer_frames = hmm_scan_parser.hmmscan_parse(folder_name + "hmmer_sequence_files/" + sequence_id + "_hmmscan_hits.fa") #only works for uniref90 as of 10/11/2021
			if (phmmer_frames == []):
				phmmer_frames_zero = "NA"
				phmmer_frames_one = "NA"
				score = "NA"
				c_e_value = "NA"
				i_e_value = "NA"
				protein_annotation_frame[-1].extend([phmmer_frames_zero, phmmer_frames_one, score, c_e_value, i_e_value])
			else:
				phmmer_frames_zero = phmmer_frames[0][0]
				phmmer_frames_one = phmmer_frames [0][1]
				score = phmmer_frames[0][2]
				c_e_value = phmmer_frames[0][3]
				i_e_value = phmmer_frames[0][4]
				protein_annotation_frame[-1].extend([phmmer_frames_zero, phmmer_frames_one, score, c_e_value, i_e_value])
		a += 1

	# write gene annotations to file:

	protein_annotation_frame_file = folder_name + "/" + frame[frame_url_index]   + "_annotation_frame.csv" # may want to add an identifier!!!!!!
	annotations_prefix = b.split("all_hits.csv")[0]
	annotations_file = annotations_prefix + "_annotations.csv"
	append_frame_url = open(annotations_file, "a")
	append_writer = csv.writer(append_frame_url)
	fcntl.flock(append_frame_url.fileno(),fcntl.LOCK_EX)
	write_frame_url = open(protein_annotation_frame_file, "w")
	spam_writer = csv.writer(write_frame_url)
	for frame in protein_annotation_frame[1:]:
		spam_writer.writerow(frame)
		append_writer.writerow(frame)
	write_frame_url.close()
	append_frame_url.close()
	shutil.rmtree(folder_name + "hmmer_sequence_files/")

	write_frame_url = open(protein_annotation_frame_file, "r")
#	fcntl.flock(write_frame_url.fileno(), fcntl.LOCK_EX)
	sql_write_frame = pandas.read_csv(write_frame_url)
	if (is_phage):
		# This will not work correctly because sets, unlike lists are not ordered!!
		header_vector = list(["phage_id","sequence_id","protein_start_site","protein_end_site","sense","orf_url","genome_length","phmmer_prediction_protein_ID","phmmer_prediction","phmmer_score","c_e-value","i_e_value","hhsearch_id","hhsearch_description","probability","E-value","score","similarity","hhblits_id","hhblits_description","hhblits_probability","hhblits_E-value","hhblits_score","hhblits_similarity","target_protein_start","target_protein_end","target_protein_sense","target_protein_order","pfam_blits_id","pfam_description","pfam_probability","pfam_E-value","pfam_score","pfam_similarity","crisprcasfinder_blits_id","crisprcasfinder_description","crisprcasfinder_score","crisprcasfinder_c_E-value","crisprcasfinder_i_E_value","type_III_signal_id","type_III_signal_description","type_III_signal_probability","type_III_signal_E_value","type_III_signal_score","type_III_signal_similarity"])
	else:
		# This will not work correctly because sets, unlike lists are not ordered!!
		header_vector = list(["genome_id","sequence_id","protein_start_site","protein_end_site","sense","orf_url","genome_length","phmmer_prediction_protein_ID","phmmer_prediction","phmmer_score","c_e-value","i_e_value","hhsearch_id","hhsearch_description","probability","E-value","score","similarity","hhblits_id","hhblits_description","hhblits_probability","hhblits_E-value","hhblits_score","hhblits_similarity","target_protein_start","target_protein_end","target_protein_sense","target_protein_order","pfam_blits_id","pfam_description","pfam_probability","pfam_E-value","pfam_score","pfam_similarity","crisprcasfinder_blits_id","crisprcasfinder_description","crisprcasfinder_score","crisprcasfinder_c_E-value","crisprcasfinder_i_E_value", "type_III_signal_id","type_III_signal_description","type_III_signal_probability","type_III_signal_E_value","type_III_signal_score","type_III_signal_similarity"])
	for header in header_vector:
		if (header not in sql_write_frame.columns):
			sql_write_frame[header] = numpy.nan # Is this meant to be for the values in the table??
	sql_write_frame.to_csv(protein_annotation_frame_file + "_annotation_frame_check.csv")
	


	# need to change code so that CRISPR-arrays are recorded as conventional rows!!
	'''
	if (is_phage):
		sql_write_frame.rename(columns={"Genome_id":"Phage_id"},inplace=True)
		cur.execute("CREATE TABLE IF NOT EXISTS PHAGE_PROTEINS(phage_id STRING, sequence_id STRING, protein_start_site STRING, protein_end_site STRING, sense STRING, orf_url STRING, genome_length STRING,  phmmer_prediction_protein_ID STRING, phmmer_prediction STRING, phmmer_score STRING, `c_e-value` STRING, i_e_value STRING, hhsearch_id STRING, hhsearch_description STRING, probability STRING, `E-value` STRING, score STRING, similarity STRING, hhblits_id STRING, hhblits_description STRING, hhblits_probability STRING, `hhblits_E-value` STRING, hhblits_score STRING, hhblits_similarity STRING,  target_protein_start STRING, target_protein_end STRING, target_protein_sense STRING, target_protein_order STRING, pfam_blits_id STRING, pfam_description STRING, pfam_probability STRING, `pfam_E-value` STRING, pfam_score STRING, pfam_similarity STRING, crisprcasfinder_blits_id STRING, crisprcasfinder_description STRING, crisprcasfinder_score STRING, `crisprcasfinder_c_E-value` STRING, `crisprcasfinder_i_E-value` STRING, padlocplus_blits_id STRING, padlocplus_description STRING, padlocplus_score STRING, `padlocplus_c_E-value` STRING,`padlocplus_i_E-value` STRING, type_III_signal_id STRING,type_III_signal_description STRING,type_III_signal_probability STRING,type_III_signal_E_value STRING,type_III_signal_score STRING,type_III_signal_similarity STRING)")
		sql_write_frame.to_sql('PHAGE_PROTEINS', con, if_exists='append', index=False)
		cur.execute("CREATE INDEX IF NOT EXISTS PHAGE_ID_INDEX ON PHAGE_PROTEINS(phage_id)")
	
	else:	
		cur.execute("CREATE TABLE IF NOT EXISTS PROTEINS(genome_id STRING, sequence_id STRING, protein_start_site STRING, protein_end_site STRING, sense STRING, orf_url STRING, genome_length STRING,  phmmer_prediction_protein_ID STRING, phmmer_prediction STRING, phmmer_score STRING, `c_e-value` STRING, i_e_value STRING, hhsearch_id STRING, hhsearch_description STRING, probability STRING, `E-value` STRING, score STRING, similarity STRING, hhblits_id STRING, hhblits_description STRING, hhblits_probability STRING, `hhblits_E-value` STRING, hhblits_score STRING, hhblits_similarity STRING,  target_protein_start STRING, target_protein_end STRING, target_protein_sense STRING, target_protein_order STRING, pfam_blits_id STRING, pfam_description STRING, pfam_probability STRING, `pfam_E-value` STRING, pfam_score STRING, pfam_similarity STRING, crisprcasfinder_blits_id STRING, crisprcasfinder_description STRING, crisprcasfinder_score STRING, `crisprcasfinder_c_E-value` STRING, `crisprcasfinder_i_E-value` STRING, padlocplus_blits_id STRING, padlocplus_description STRING, padlocplus_score STRING, `padlocplus_c_E-value` STRING,`padlocplus_i_E-value` STRING, type_III_signal_id STRING,type_III_signal_description STRING,type_III_signal_probability STRING,type_III_signal_E_value STRING,type_III_signal_score STRING,type_III_signal_similarity STRING)")
		sql_write_frame.to_sql('PROTEINS', con, if_exists='append', index=False)
		cur.execute("CREATE INDEX IF NOT EXISTS GENOME_ID_INDEX ON PROTEINS(genome_id)")
	'''
	# unused seed sequence code
	'''
	if (target_pos_switch == 1):
		# add option to retrieve this information from corresponding csv file as well as sql table!!!!!!!!!!!!!!!!!!
		# need to load all_hits_table + spacer_hit_map table. Dictionalise on Genome_id and Phage id
		for protein in protein_annotation_frame: # should be the same
			if (is_phage == False):
				my_row = cur.execute("SELECT match_start, match_end FROM GENOMES WHERE Genome_id=?",(protein[0],))
				startend = my_row.fetchone()
			else:
				my_row = cur.execute("SELECT Mapped_start_site, Mapped_end_site FROM SPACER_HITMAP WHERE Phage_id=?",(protein[0],))
				startend = my_row.fetchone()

		# check that the mapped start/end sites are in ascending order
			if (startend[0] > startend[1]):
				startend_cpy = copy.deepcopy(startend[0])
				startend[0] = startend[1]
				startend[1] = startend_cpy
				if (is_phage):
					break
		# need to align these coordinates with the phage coordinates
			if (is_phage == False):
				protein_start = int(protein[2])
				protein_end = int(protein[3])
				protein_sense = int(protein[4])
				if (protein_start < startend[0] < protein_end or protein_start < startend[1] < protein_end):
					startend[0] = protein_start
					startend[1] = protein_end
					break
	'''				
	con.commit()
	write_frame_url.close()
	subprocess.run(["rm -r " + folder_name],shell=True)
	subprocess.run(["rm *.ps"],shell=True)
	return 0	



