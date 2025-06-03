from Bio import SeqIO
import sys
import csv


# This script was not used. ncbi_esearch_final_annotation_retriever5.py was used instead.
# INPUT: 1. table/list containing NCBI sequence identifers
host_url = open(sys.argv[1],"r")
host_table = csv.reader(host_url)

host_ids = {}

for row in host_table:
	my_id = row[0].split("::")[0]
	if my_id not in host_ids:
		host_ids[my_id] = my_id

host_url.close()

host_ids = list(host_ids.keys())
# INPUT 2. master genbank table from NCBI FTP repository with additional entries added by entrez queries.
# i.e.master_genbank_table.csv
ncbi_url = open(sys.argv[2],"r")
ncbi_table = csv.reader(ncbi_url)

ncbi_dict = {}
for row in ncbi_table:
	ncbi_dict[row[0].split("::")[0]] = row

ncbi_url.close()
print(sys.argv[1])
print(len(host_ids))
in_count = 0
ret_out = open(sys.argv[3],"w")
spamwriter = csv.writer(ret_out)
# OUTPUT: Table containing matching annotations from NCBI database for the subset of sequences specified in the input.
unmatched_out = open(sys.argv[3] + "_unmatched.csv","w")
sp_out = csv.writer(unmatched_out)
out_count = 0
for name in host_ids:
	my_name = name
	if (my_name in ncbi_dict):
		in_count += 1 
		spamwriter.writerow(ncbi_dict[my_name])
	else:
		sp_out.writerow([name])
		out_count += 1
ret_out.close()
unmatched_out.close()		
print(in_count, out_count)	


