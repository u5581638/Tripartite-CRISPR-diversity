# properties_species_annotation_functions.py
import math
import re
# Script to hold all helper functions used for plasmid/virsort/species annotation

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
		if (row[27] == "NCLDV"):
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
	if (denominator > 0):		
		genome_linear = genome_linear / denominator
		genome_circular = genome_circular / denominator
	else:
		genome_linear = 0
		genome_circular = 0	
	
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
		print(entry)
		print(index)
		print("good")
		my_entry = entry[index].split(" ")
		print(my_entry)
		# This syntax should be simplified if possible!!
		if (len(entry) > 1 and entry[index+1] != '' and my_entry[-1] == "metagenome" and len(my_entry) > 1):
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
		gold_id = gene_id.split("_")[0]
		gene_id = re.split('[0-9]+', gene_id) [0]
		if (gold_id in gold_anno_dict):
			matching_gold_entries.append(gold_anno_dict[gold_id])

		elif (gene_id in ncbi_dict):
			matching_ncbi_entries.append(ncbi_dict[gene_id])
		else:
			print("missing")
			print(gene_id)	
			print(gold_id)

	gold_sample_metagenome = 0
	gold_sample_bacteria = 0
	gold_sample_archaea = 0
	gold_sample_phage = 0
	gold_sample_unassigned = 0
	gold_species_no = 0

	if (len(matching_gold_entries) > 0):
	#	print(matching_gold_entries)
		sample_list = species_count(gold_species_clean(matching_gold_entries, 26),26) # this is the problematic line!!
		sample_list.sort(key=sp_sort,reverse=True)
		highest_ele = sample_list[0][-1]
		highest_sample_gold_abundance = highest_ele / len(matching_gold_entries)
	#	print("sample_list:")
	#	print(sample_list)
		highest_sample_gold_diversity = shannon_index(sample_list)

		species_list = species_count(gold_species_clean(matching_gold_entries, 16),16)
	
		species_list.sort(key=sp_sort,reverse=True)
		highest_ele = species_list[0][-1]

		highest_species_gold_abundance = highest_ele / len(matching_gold_entries)
		highest_species_gold_diversity = shannon_index(species_list)
	else:
		highest_sample_gold_abundance = 0
		highest_sample_gold_diversity = 0
		highest_species_gold_abundance = 0
		highest_species_gold_diversity = 0	
		print("Bugger1")
		print(cluster)

	if (len(matching_ncbi_entries) > 0):
		# This might not work/will work differently!!
		sample_ncbi_list = species_count(ncbi_species_clean(matching_ncbi_entries,25),25)
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
		print("bugger2")
		print(cluster)

	return (highest_sample_gold_abundance,highest_sample_gold_diversity,highest_sample_ncbi_abundance,highest_sample_ncbi_diversity,highest_species_gold_abundance,highest_species_gold_diversity,highest_species_ncbi_abundance,highest_species_ncbi_diversity)
		 #need to figure what best to return. Not clear what to return from species info directly
def species_survey (cluster,gold_anno_dict,ncbi_dict):
	# need to split both species and sample into multiple categories
	# GOLD database == row 8
	# ncbi databse == row 26
	matching_gold_entries = []
	matching_ncbi_entries = []

	for gene_id in cluster:
	#	gene_id = gene_id.split(":")[0]
		if (gene_id in gold_anno_dict):
			matching_gold_entries.append(gold_anno_dict[gene_id])
		
		# need to change the way this is formulated.
		
		if (gene_id in ncbi_dict):
			matching_ncbi_entries.append(ncbi_dict[gene_id])

	if (len(matching_gold_entries) > 0):
	#	print(matching_gold_entries)
		sample_list = species_count(gold_species_clean(matching_gold_entries, 28),28) # this is the problematic line!!
		sample_list.sort(key=sp_sort,reverse=True)
		species_list = species_count(gold_species_clean(matching_gold_entries, 15),15)	
		species_list.sort(key=sp_sort,reverse=True)

	if (len(matching_ncbi_entries) > 0):
		# This might not work/will work differently!!
		sample_ncbi_list = species_count(ncbi_species_clean(matching_ncbi_entries,25),25)
		sample_ncbi_list.sort(key=sp_sort,reverse=True)
		species_ncbi_list = species_count(ncbi_species_clean(matching_ncbi_entries,8),8)
		species_ncbi_list.sort(key=sp_sort,reverse=True)
	
	if (len(matching_gold_entries) > 0 and len(matching_ncbi_entries) > 0):
		final_list = species_list + sample_list + species_ncbi_list + sample_ncbi_list 
		final_list.sort(key=sp_sort,reverse=True)
		total_size = 0
		for ele in final_list:
			if (ele[0] != ""):
				total_size += int(ele[1])

		return (total_size,final_list)
	elif (len(matching_gold_entries) > 0):
		final_list = species_list + sample_list
		final_list.sort(key=sp_sort,reverse=True)
		total_size = 0
		for ele in final_list:
			if (ele[0] != ""):
				print(ele)
				total_size += int(ele[1])
		return (total_size,final_list)
	elif (len(matching_ncbi_entries) < 0):
		final_list = species_ncbi_list + sample_ncbi_list
		final_list.sort(key=sp_sort,reverse=True)
		total_size = 0
		for ele in final_list:
			if (ele[0] != ""):
				total_size += int(ele[1])
		return (total_size,final_list)
	else:
		return (0,[])

def matrix_generation(component_list):
	out_matrix  = []
	header = set()
	for component in component_list:
		cluster = component[3]
		for pair in cluster:
			if pair[0] not in header:
				header.add(pair[0])
	header = list(header)	
	out_matrix.append(["Component_no","proportion_detected"] + header)		
	for component in component_list:
		comp_dict = {}
		cluster = component[3]
		for pair in cluster:
			if (pair[0] not in comp_dict):
				comp_dict[pair[0]] = pair[1]
			else:
				print("comp1_error!!")	
		outrow = []
		outrow.append(component[0])
		print(component)
		outrow.append(float(component[1]) / float(component[2]))
		for head in header:
			if (head in comp_dict):
				outrow.append(float(comp_dict[head]) / float(component[2]))
			else:
				outrow.append("0")
		out_matrix.append(outrow)
	return out_matrix





	# first iterate through all components and sort in descending order.
	# list names as a set. 
	# for each component assign according to where the label is in the list elements
	# consider adding a special row called total protein to denote the total fraction detected.