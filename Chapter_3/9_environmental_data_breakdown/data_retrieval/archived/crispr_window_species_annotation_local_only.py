# crispr_window_species_annotation

from Bio import SeqIO
import sys
import csv

# INPUT: 1. DNA extracted within a 40kb window of the CRISPR array ("spliced_debug_corrected_cleaned_rerun_20kb_local_windows_combined.fasta")
#		 2. merged gold annotation table ("true_gold_compile3_merged.csv")
# OUTPUT: table of JGI-matched metadata entries for crispr-windows

sequences = SeqIO.parse(sys.argv[1],"fasta")

with open(sys.argv[2], "r") as csvfile:
	hit_table = csv.reader(csvfile)
	hit_dict = {}
	for row in hit_table:
		if (row[0] not in hit_dict):
			hit_dict[row[0]] = row
		else:
			print("Error!!")


ret_out = open(sys.argv[3],"w")
spamwriter = csv.writer(ret_out)
spamwriter.writerow(["Sequence_writer","AP_id","AP_NAME","AP_GENBANK","AP_SRA","AP_GOLD_PROJECT_ID","PROJECT_NAME","SEQUENCING_STRATEGY","STUDY_ID","ORGANISM_ID","BIOSAMPLE_ID","ORGANISM_NAME","PHYLUM","CLASS","ORDER","FAMILY","GENUS","SPECIES","GRAM","BIOSAMPLE_NAME","NAME","COLLECTION_SITE","LOCATION","BIOSAMPLE_ECOSYSTEM","ECOSYSTEM_TYPE","ECOSYSTEM_SUBTYPE","ECOSYSTEM"])
missed_list = open(sys.argv[3] + "_missed.csv","w")
missed_writer = csv.writer(missed_list)
for sequ in sequences:
	sequ_id = sequ.id.split("::") [0]
	sequ_id = sequ_id.split("_") [0]
#	print(sequ_id)
#	print(hit_dict.keys())
#	exit()
	if sequ_id in hit_dict:
		spamwriter.writerow([sequ.id] + hit_dict[sequ_id])
	else:
		missed_writer.writerow([sequ.id])




ret_out.close()
missed_list.close()