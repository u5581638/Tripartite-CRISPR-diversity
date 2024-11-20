# program to generate a matrix to be used as the input to an R based heatmap generation function to summarise and represent substructure differences between different input metagenomic sequences
# Inputs:

# need to input the network. Both host and phage

# Input as host-phage table precursor.
# Need to compensate for mmseqs clustering in cases where this is required


# Host/Phage
# 1. Virsorter
# 1. a) Phage identity?
# 1. b) Phage shape (linear/circular)
# 1. c) Sequence type (host/prophage/free phage/plasmid)
# 2. Plasme
# 2. a) Phage shape (linear/circular)
# 2. b) Sequence type (host/prophage/free phage/plasmid)
# 3. JGI and NCBI annotations
# 3. a) sample info (if applicable - cultured samples have no explicit info)
# 3. b) species info
# 4. annotated protein composition (host only)
# 4. a) use PFAM
# 5. PAM identity by clade (host only)

# heatmaps can only represent these properties in binary ways.

# Outline:
# 1. Phage species highest abundance. (define this: highest occurrance of the most common species / total occurance in the associated cluster )
# 2. Phage species high abundance / diversity (define this: highest occurrance of the most common species / number of different species annotations)
# 3. Phage circular
# 4. Phage linear
# 5. host
# 6. prophage
# 7. free phage
# 8. plasmid
# 9.  Phage sample highest abundance. (define this: highest occurrance of the most common species / total occurance in the associated cluster )
# 10. Phage sample high abundance / diversity (define this: highest occurrance of the most common species / number of different species annotations)


import sys
import csv
import math

# 1. import host-phage interaction table



def species_count (species_rows,index):
#	print("rows")
#	print(species_rows)
	ret_dict = {}
	for row in species_rows:
		if (row not in ret_dict):
			ret_dict[row] = [row,1]
		else:
			ret_dict[row][1] += 1	
	
	return list(ret_dict.values())	

def virsort_proportion (cluster,virsort_dict):
	matching_entries = []
	#print(len(cluster))
	for gene_id  in cluster:
		if gene_id in virsort_dict:
			matching_entries.append(virsort_dict[gene_id])

	genome_linear = 0
	genome_circular = 0
	cell = 0
	prophage = 0 # unsure if this is detectable
	phage = 0
	dsDNAphage = 0
	ssDNAphage = 0
	rna_phage = 0
	lavidaviridae = 0
	ncldv = 0
	denominator = len(cluster)
	for row in matching_entries:
	#	print("my_row:")
		print(row[27],row[28])
		if (row[28] == "linear"):
			genome_linear += 1
		if (row[28] == "circular"):
			genome_circular += 1
		# possible that no entry is given if no virus is detected!
		if (row[27] == "dsDNAphage"):
			dsDNAphage += 1
		if (row[27] == "ssDNA"):
			ssDNAphage += 1
		if (row[27] == "RNA"):
			rna_phage += 1
		if (row[27] == "lavidaviridae"):
			lavidaviridae += 1
		if (row[27] == "ncldv"):
			ncldv += 1

	if (denominator > 0):
	#	print("denominator:")
	#	print(denominator)

		genome_linear = genome_linear / denominator
		genome_circular = genome_circular / denominator
		cell = 1 - ((dsDNAphage + ssDNAphage + rna_phage + lavidaviridae + ncldv) / denominator)
		dsDNAphage = dsDNAphage / denominator
		ssDNAphage = ssDNAphage / denominator	
		rna_phage = rna_phage / denominator
		lavidaviridae = lavidaviridae / denominator
		ncldv = ncldv / denominator
	else:
		genome_linear = 0
		genome_circular = 0
		dsDNAphage = 0
		ssDNAphage = 0
		cell = 0
		rna_phage = 0
		lavidaviridae = 0
		ncldv = 0

	# shape = row 28 
	# genome type = row 29		
	return (genome_linear, genome_circular, cell, dsDNAphage, ssDNAphage, rna_phage, lavidaviridae, ncldv)

def plasme_proportion (cluster,plasme_dict):
	matching_entries = []
	for gene_id in  cluster:
		if (gene_id in plasme_dict):
			matching_entries.append(plasme_dict[gene_id])
	# finish writing code here!!
	genome_linear = 0
	genome_circular = 0
	denominator = len(matching_entries)
	for row in matching_entries:
		if (row[2] == "None"):
			genome_linear += 1
		else: # if (row[2] == "circular")
			genome_circular += 1
	genome_linear = genome_linear / denominator
	genome_circular = genome_circular / denominator
	
	return (genome_linear,genome_circular)

def sp_sort(x):
	return x[1]

# This metric is postive when the species are diverse and negative when there are more members of a small number of species and few species types.
def shannon_index(species): # N, list(n)
	individual_count = 0
#	print("species:")
#	print(len(species))
	total_count = 0
	for i in species:
		total_count += i[1]
		pi = i[1]
	#	print(pi)
		individual_count += (pi * math.log(pi))
#	print(individual_count)
	if (total_count > 0):
		h_index = (total_count * math.log(total_count) - individual_count) / (total_count)
	else:
		h_index = 1	
	if (len(species) == 1):	
		return 1 # 0 / 0 = 1 in this context!!
	else:	
		return h_index / math.log(len(species))

def ncbi_species_clean(matching_entries_list, index):
	ret_list = []
	for entry in matching_entries_list:
		if (len(entry) >= index):
			my_entry = entry[index].split(" ")
			if (my_entry[-1] == "metagenome" and len(my_entry) > 1):
				my_entry = my_entry[0] + " " + my_entry[1]
			else:
				my_entry = my_entry[0]
		#	entry[8] = my_entry
			ret_list.append(my_entry)
		else:
			continue		

	return ret_list
# attempt to convert the first 1-2 words of the species id into a key based on whether data is metagenome only or contains species information
def gold_species_clean(matching_entries_list,index):	
	ret_list = []
	for entry in matching_entries_list:

		my_entry = entry[index].split(" ")
		# This syntax should be simplified if possible!!
		if (len(entry) > 1 and entry[index+1] != '' and entry[index+1][-1] == "metagenome" and len(my_entry) > 1):
			my_entry = my_entry[0] + " " + my_entry[1]
		else:
			my_entry = my_entry[0]

		ret_list.append(my_entry)
	return ret_list	

def species_annotation (cluster,gold_anno_dict,ncbi_dict):
	# need to split both species and sample into multiple categories
	# GOLD database == row 8
	# ncbi databse == row 26
	matching_gold_entries = []
	matching_ncbi_entries = []

	for gene_id in cluster:
		gene_id = gene_id.split(":")[0]
		if (gene_id in gold_anno_dict):
			matching_gold_entries.append(gold_anno_dict[gene_id])
		
		# need to change the way this is formulated.
		
		if (gene_id in ncbi_dict):
			matching_ncbi_entries.append(ncbi_dict[gene_id])

	gold_sample_metagenome = 0
	gold_sample_bacteria = 0
	gold_sample_archaea = 0
	gold_sample_phage = 0
	gold_sample_unassigned = 0
	gold_species_no = 0

	if (len(matching_gold_entries) > 0):
	#	print(matching_gold_entries)
		sample_list = species_count(gold_species_clean(matching_gold_entries, 28),28) # this is the problematic line!!
		sample_list.sort(key=sp_sort,reverse=True)
		highest_ele = sample_list[0][-1]
		highest_sample_gold_abundance = highest_ele / len(matching_gold_entries)
	#	print("sample_list:")
	#	print(sample_list)
		highest_sample_gold_diversity = shannon_index(sample_list)

		species_list = species_count(gold_species_clean(matching_gold_entries, 15),15)
	
		species_list.sort(key=sp_sort,reverse=True)
		highest_ele = species_list[0][-1]

		highest_species_gold_abundance = highest_ele / len(matching_gold_entries)
		highest_species_gold_diversity = shannon_index(species_list)
	else:
		highest_sample_gold_abundance = 0
		highest_sample_gold_diversity = 0
		highest_species_gold_abundance = 0
		highest_species_gold_diversity = 0	

	if (len(matching_ncbi_entries) > 0):
		# This might not work/will work differently!!
		my_clean = ncbi_species_clean(matching_ncbi_entries,25)
		sample_ncbi_list = species_count(my_clean,25)
		sample_ncbi_list.sort(key=sp_sort,reverse=True)
		highest_ncbi_ele = sample_ncbi_list[0][-1]
		highest_sample_ncbi_abundance = highest_ncbi_ele / len(matching_ncbi_entries)
		highest_sample_ncbi_diversity = shannon_index(sample_ncbi_list)

		species_ncbi_list = species_count(ncbi_species_clean(matching_ncbi_entries,8),8)
		species_ncbi_list.sort(key=sp_sort,reverse=True)
		highest_ncbi_ele = species_ncbi_list[0][-1]
		highest_species_ncbi_abundance = highest_ncbi_ele / len(matching_ncbi_entries)
		highest_species_ncbi_diversity = shannon_index(species_ncbi_list)
	else:
		highest_sample_ncbi_abundance = 0
		highest_sample_ncbi_diversity = 0
		highest_species_ncbi_abundance = 0
		highest_species_ncbi_diversity = 0

	return (highest_sample_gold_abundance,highest_sample_gold_diversity,highest_sample_ncbi_abundance,highest_sample_ncbi_diversity,highest_species_gold_abundance,highest_species_gold_diversity,highest_species_ncbi_abundance,highest_species_ncbi_diversity)
	 #need to figure what best to return. Not clear what to return from species info directly

# INPUTS:
# 1.:	 .clusters file generated via vCONTACT2
with open(sys.argv[1],"r") as csvfile1:
	host_cluster_table = list(csv.reader(csvfile1)) 
	# may need to add a table to map genome_ids to cluster_ids from vContact2


# with open(sys.argv[2],"r") as csvfile1b:
#	phage_cluster_table = list(csv.reader(csvfile1b))
# host files
# 2.:	"final_table.csv" in virsorter output folder (from annotation pipeline)
with open(sys.argv[2],"r") as csvfile2:
	host_virsort_table = list(csv.reader(csvfile2))

virsort_dict = {}

for row in host_virsort_table:
	if (row[1] not in virsort_dict):
		virsort_dict[row[1]] = row
			
# 3.:
with open(sys.argv[3],"r") as csvfile3:
	host_plasme_table = list(csv.reader(csvfile3))

plasme_dict = {}

for row in host_plasme_table:
	if (row[0] not in plasme_dict):
		plasme_dict[row[0]] = row	




with open(sys.argv[4],"r") as csvfile4:
	host_gold_annotation_table = list(csv.reader(csvfile4)) # need to a script to create this file. Must be specific to a CRISPR-subtype.

gold_dict = {}
for row in host_gold_annotation_table[1:]:
	if row[0] not in gold_dict:
		gold_dict[row[0]] = row

# Theres's a problem with the script that generates this input file.
with open(sys.argv[5],"r") as csvfile5:
	host_ncbi_annotation_table = list(csv.reader(csvfile5)) # need to a script to create this file. Must be specific to a CRISPR-subtype.


ncbi_dict = {}
for row in host_ncbi_annotation_table:
	if row[0] not in ncbi_dict:
		ncbi_dict[row[0]] = row	

# need a seperate heatmap/annotation structure to show protein annotations by cluster
# need a seperate heatmap/annotation structure to show PAMs by cluster.

with open(sys.argv[6],"r") as csvfile:
	component_table = list(csv.reader(csvfile))

component_dict = {}
for comp in component_table[1:]:
	if (comp[1] not in component_dict):
		component_dict[comp[1]] = comp[2]
	else:
		print("error") 

output_matrix = []
host_phage_prefix = "p"
# 1. need to compute the occurrance of each trait
# 2. need to deal with each host sequence by cluster

components = {}

cluster_number = 1	
# This loop will need to be split into two
for row in host_cluster_table[1:]:
	# is this a component
	if (host_phage_prefix + str(cluster_number) in component_dict):
	 	component_number = component_dict[host_phage_prefix + str(cluster_number)]
#	 	print("Yay!!")
	 	if (component_number not in components):
	 		components[component_number] = row[7].split(" ")
	 	else:
	 		my_rows = row[7].split(" ")
	 		components[component_number].extend(my_rows)
	else:
		print("error, all components should be in here")
	cluster_number += 1	 
cluster_number = 1	
for component in components:
	# need to assign arbitary cluster numbers
	
	 # want to aggregate cluster members from all components into a 
	 cluster_members = components[component]
	 virsort_properties = virsort_proportion(cluster_members,virsort_dict)
	 plasme_properties = plasme_proportion(cluster_members,plasme_dict)
	 species_properties = species_annotation(cluster_members,gold_dict,ncbi_dict)
	 # Fill this in after running the main functions tomorrow!!
	 # The component number should be the key
	 output_matrix.append(["component_" + str(component)] + list(virsort_properties) + list(plasme_properties) + list(species_properties))

my_out = open(sys.argv[7],"w")
spamwriter = csv.writer(my_out)
spamwriter.writerow(["cluster_no","virsort_linear","virsort_circular","virsort_cell","virsort_dsDNAphage","virsort_ssDNAphage","virsort_RNAphage","virsort_lavidaviridae","virsort_ncldv" ,"plasme_linear","plasme_circular", "gold_sample_abundance","gold_sample_diversity","highest_sample_ncbi_abundance","highest_sample_ncbi_diversity","highest_genus_gold_abundance","highest_genus_gold_diversity","highest_genus_ncbi_abundance","highest_genus_ncbi_diversity"])
for row in output_matrix:
	spamwriter.writerow(row)
my_out.close()




