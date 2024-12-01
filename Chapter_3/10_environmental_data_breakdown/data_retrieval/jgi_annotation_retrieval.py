from Bio import SeqIO
import sys
import csv

# script to retrieve sequence metadata for JGI-derived sequences using a merged metadata table downloaded from the GOLD database website

# INPUT: table containing JGI sequence identifers
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
jgi_url = open(sys.argv[2],"r")
jgi_table = csv.reader(jgi_url)

jgi_dict = {}
for row in jgi_table:
	jgi_dict[row[12]] = row

jgi_url.close()
print(sys.argv[1])
print(len(host_ids))
in_count = 0
ret_out = open(sys.argv[3],"w")
spamwriter = csv.writer(ret_out)
unmatched_out = open(sys.argv[3] + "_unmatched.csv","w")
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


