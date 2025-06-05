# host-phage vCONTACT clade component annotation
# script to annotate the components of vCONTACT2 generated networks by number of species.
# Should work with p/h cluster numbers so that the chart can be compared to the bipartite phage-host interaction diagrams.
# Steps:
# 1. Load host/phage vCONTACT2 .clusters + .ntw file + other annotations.
# 2. Find the number of components in the .ntw file
# 3. Find the annotation proportion of eac component. This may require developing new functions to retrieve the species.
# 4. Use species function to retrieve annotations
# 5. Need to label the graph diagram with component numbers. <- this will also be the basis for overlaying other stuff.
import sys
import csv
import math
from properties_species_annotation_functions import virsort_proportion
from properties_species_annotation_functions import plasme_proportion
from properties_species_annotation_functions import species_survey
from properties_species_annotation_functions import matrix_generation

# INPUT: 
# 1. host/phage c1.clusters file in csv format (from vConTACT output)
# make this the network file. Extract directly from the network.
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
			
# 3.: table of genome shape predictions returned by PLasME.
with open(sys.argv[3],"r") as csvfile3:
	host_plasme_table = list(csv.reader(csvfile3))

plasme_dict = {}

for row in host_plasme_table:
	if (row[0] not in plasme_dict):
		plasme_dict[row[0]] = row	

# 4. Table of JGI metadata annotations (from GOLD) i.e. 
with open(sys.argv[4],"r") as csvfile4:
	host_gold_annotation_table = list(csv.reader(csvfile4)) # need to a script to create this file. Must be specific to a CRISPR-subtype.

gold_dict = {}
for row in host_gold_annotation_table[1:]:
	if row[0] not in gold_dict:
		gold_dict[row[0]] = row

# 5. Table of ncbi metadata annotations
with open(sys.argv[5],"r") as csvfile5:
	host_ncbi_annotation_table = list(csv.reader(csvfile5)) # need to a script to create this file. Must be specific to a CRISPR-subtype.


ncbi_dict = {}
for row in host_ncbi_annotation_table:
	if row[0] not in ncbi_dict:
		ncbi_dict[row[0]] = row	

# 6. table of sequence identifiers and partitioned cluster numbers.
with open(sys.argv[6],"r") as csvfile:
	component_table = list(csv.reader(csvfile))

# OUTPUT: Conservation matrix showing the proportion of linear/circular+host/phage(RNA/dsDNA/ssDNA/jumbo phage) in each partitioned cluster
# i.e cas13b_virsort_phage_matrix.csv

component_dict = {}
for comp in component_table[1:]:
	if (comp[2] not in component_dict):
		component_dict[comp[2]] = [comp[1]]
	else:
		component_dict[comp[2]].append(comp[1])

output_matrix = []
host_phage_prefix = "p"
# 1. need to compute the occurrance of each trait
# 2. need to deal with each host sequence by cluster

components = component_dict

component_list = []
for component in components:
	# need to assign arbitary cluster numbers
	
	 # want to aggregate cluster members from all components into a 
	 cluster_members = components[component]
	 virsort_properties = virsort_proportion(cluster_members,virsort_dict)
	 plasme_properties = plasme_proportion(cluster_members,plasme_dict)
	# species_properties = species_annotation(cluster_members,gold_dict,ncbi_dict)
	# may wish to convert to proportions of the component or component detected!
	# species_properties = species_survey(cluster_members,gold_dict,ncbi_dict)
	# cluster_list = species_properties[1]
	 # combined size of matching entries
	# total_size = species_properties[0]
	 # combined size of sequences in component
	# component_total_size = len(cluster_members)
	 # Fill this in after running the main functions tomorrow!!
	 component_list.append(["component_" + str(component)] + list(virsort_properties) + list(plasme_properties)) # may wish to add virsort + plasme properties
	 # The component number should be the key
#	 print(species_properties)
#	 output_matrix.append(["component_" + str(component)] + list(virsort_properties) + list(plasme_properties) + list(species_properties))

# need to reconcile these properties into a single set of matrix columns
output_matrix =component_list

my_out = open(sys.argv[7],"w")
spamwriter = csv.writer(my_out)
# Should not need this header line anymore
spamwriter.writerow(["Component_no","virsort_linear","virsort_circular","virsort_cell","virsort_dsDNAphage","virsort_ssDNAphage","virsort_RNAphage","virsort_lavidaviridae","virsort_ncldv" ,"plasme_linear","plasme_circular"])

for row in output_matrix:
	spamwriter.writerow(row)
my_out.close()

