from Bio import SeqIO
import sys
import csv

# script to retrieve sequence metadata for JGI-derived sequences using a merged metadata table downloaded from the GOLD database website

# SHELL: find labelled_genomes -name "*.fasta" | xargs -n 1 -I {} -P 1 python3 jgi_annotation_retrieval.py {} true_gold_compile3_merged.csv {}_gold_annotations.csv
# INPUT: table containing JGI sequence identifers (Ga numbers)
host_url = open(sys.argv[1],"r")
host_table = csv.reader(host_url)

host_ids = {}

for row in host_table:
	my_id = row[0].split("::")[0]
	if my_id not in host_ids:
		host_ids[my_id] = my_id

host_url.close()

host_ids = list(host_ids.keys())
# INPUT: GOLD table (merged)
# i.e. true_gold_compile3_merged.csv
jgi_url = open(sys.argv[2],"r")
jgi_table = csv.reader(jgi_url)

jgi_dict = {}
for row in jgi_table:
	jgi_dict[row[12]] = row

jgi_url.close()
print(sys.argv[1])
print(len(host_ids))
in_count = 0

# OUTPUT: Table containing matching annotations from JGI database for the subset of sequences specified in the input.
ret_out = open(sys.argv[3],"w")
spamwriter = csv.writer(ret_out)
# changed name
unmatched_out = open(sys.argv[3] + "_missed.csv","w")
sp_out = csv.writer(unmatched_out)
out_count = 0
for name in host_ids:
	my_name = name.split("_") [0]
	if (my_name in jgi_dict):
		in_count += 1 
		spamwriter.writerow(jgi_dict[my_name])
	else:
		sp_out.writerow([name])
		out_count += 1
ret_out.close()
unmatched_out.close()		
print(in_count, out_count)	


