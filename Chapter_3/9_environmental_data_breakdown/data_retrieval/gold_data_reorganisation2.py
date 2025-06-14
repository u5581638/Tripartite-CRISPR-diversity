# gold data reorganisation
# script to reorganise the GOLD table into a single table indexed by analysis id showing:
# 1. biosample information
# 2. species information
# 3. accessions (if applicable - SRA is most desired)

import csv
import sys

with open(sys.argv[1],"r") as csvfile:
	biosample_table = list(csv.reader(csvfile))
	biosample_dict = {}
	for sample in biosample_table:
		sample_id = sample[0]
		if sample_id not in biosample_dict:
			biosample_dict[sample_id] = sample

with open(sys.argv[2],"r") as csvfile2:
	organism_table = list(csv.reader(csvfile2))
	organism_dict = {}
	for organism in organism_table:
		organism_id = organism[0]
		if (organism_id not in organism_dict):
			organism_dict[organism_id] = organism

with open(sys.argv[3],"r") as csvfile3:
	sequence_project_table = list(csv.reader(csvfile3))
	sp_dict = {}
	for sequence_project in sequence_project_table:
		sp_id = sequence_project[0]
		if (sp_id not in sp_dict):
			sp_dict[sp_id] = sequence_project

with open(sys.argv[4],"r") as csvfile4:
	analysis_table = list(csv.reader(csvfile4))

out_table = open(sys.argv[5],"w")
spamwriter = csv.writer(out_table)
spamwriter.writerow(["AP_id","AP_NAME","AP_GENBANK","AP_SRA","AP_GOLD_PROJECT_ID","PROJECT_NAME","SEQUENCING_STRATEGY","STUDY_ID","ORGANISM_ID","BIOSAMPLE_ID","ORGANISM_NAME","PHYLUM","CLASS","ORDER","FAMILY","GENUS","SPECIES","GRAM","ORGANISM_HABITAT","ORGANISM_COLLECTION_SITE","ORGANISM_ECOSYSTEM","ORGANISM_ECOSYSTEM_CATEGORY","ORGANISM_ECOSYSTEM_TYPE","ORGANISM_ECOSYSTEM_SUBTYPE","ORGANISM_SPECIFIC_ECOSYSTEM","ORGANISM_CULTURED","BIOSAMPLE_NAME","NAME","COLLECTION_SITE","LOCATION","BIOSAMPLE_ECOSYSTEM","ECOSYSTEM_TYPE","ECOSYSTEM_SUBTYPE","ECOSYSTEM"])
for ap in analysis_table[1:]:
	out_row = []
	out_row.append(ap[0])
	out_row.append(ap[1])
	out_row.append(ap[6])
	out_row.append(ap[7])
	out_row.append(ap[10])
	if (ap[10] in sp_dict):
		out_row.append(sp_dict[ap[10]][1])
		out_row.append(sp_dict[ap[10]][2])
		out_row.append(sp_dict[ap[10]][17])
		out_row.append(sp_dict[ap[10]][18])
		out_row.append(sp_dict[ap[10]][19])
		if (sp_dict[ap[10]][18] in organism_dict):
			out_row.append(organism_dict[sp_dict[ap[10]][18]][1])
			out_row.append(organism_dict[sp_dict[ap[10]][18]][5])
			out_row.append(organism_dict[sp_dict[ap[10]][18]][6])
			out_row.append(organism_dict[sp_dict[ap[10]][18]][7])
			out_row.append(organism_dict[sp_dict[ap[10]][18]][8])
			out_row.append(organism_dict[sp_dict[ap[10]][18]][9])
			out_row.append(organism_dict[sp_dict[ap[10]][18]][10])
			out_row.append(organism_dict[sp_dict[ap[10]][18]][15])
			out_row.append(organism_dict[sp_dict[ap[10]][18]][17])
			out_row.append(organism_dict[sp_dict[ap[10]][18]][18])
			out_row.append(organism_dict[sp_dict[ap[10]][18]][20])
			out_row.append(organism_dict[sp_dict[ap[10]][18]][21])
			out_row.append(organism_dict[sp_dict[ap[10]][18]][22])
			out_row.append(organism_dict[sp_dict[ap[10]][18]][23])
			out_row.append(organism_dict[sp_dict[ap[10]][18]][24])
			out_row.append(organism_dict[sp_dict[ap[10]][18]][43])
		else:
			out_row.append("")
			out_row.append("")
			out_row.append("")
			out_row.append("")
			out_row.append("")
			out_row.append("")
			out_row.append("")
			out_row.append("")
			out_row.append("")
			out_row.append("")
			out_row.append("")
			out_row.append("")
			out_row.append("")
			out_row.append("")
			out_row.append("")
			# add empty entries to outrow	
		if (sp_dict[ap[10]][19] in biosample_dict):
			out_row.append(biosample_dict[sp_dict[ap[10]][19]][1])
			out_row.append(biosample_dict[sp_dict[ap[10]][19]][3])
			out_row.append(biosample_dict[sp_dict[ap[10]][19]][4])
			out_row.append(biosample_dict[sp_dict[ap[10]][19]][5])
			out_row.append(biosample_dict[sp_dict[ap[10]][19]][9])
			out_row.append(biosample_dict[sp_dict[ap[10]][19]][10])
			out_row.append(biosample_dict[sp_dict[ap[10]][19]][11])
			out_row.append(biosample_dict[sp_dict[ap[10]][19]][12])
			out_row.append(biosample_dict[sp_dict[ap[10]][19]][13])
		else:	
			out_row.append("")
			out_row.append("")
			out_row.append("")
			out_row.append("")
			out_row.append("")
			out_row.append("")
			out_row.append("")
			out_row.append("")
			out_row.append("")
	else:
		out_row.append("")
		out_row.append("")
		out_row.append("")
		out_row.append("")
		out_row.append("")
		out_row.append("")
		out_row.append("")
		out_row.append("")
		out_row.append("")
		out_row.append("")
		out_row.append("")
		out_row.append("")
		out_row.append("")
		out_row.append("")
		out_row.append("")
		out_row.append("")
		out_row.append("")
		out_row.append("")
		out_row.append("")
		out_row.append("")
		out_row.append("")
		out_row.append("")
		out_row.append("")
		out_row.append("")
		out_row.append("")
		out_row.append("")
		out_row.append("")
		out_row.append("")
		out_row.append("")
		out_row.append("")
	spamwriter.writerow(out_row)


	# append empty rows for sequence project + organism + biosample

out_table.close()