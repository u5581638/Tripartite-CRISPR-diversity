from Bio import SeqIO
import sys
import csv

# do two things
# 1. Retrieve all matching JGI scaffolds from GOLD records.
# 2. Compute matched/unmatched proportion

# INPUT: 1. list of host genome ids from each subtype.
#		 2. merged_table containing JGI sequence identifers (Ga numbers). This must be: true_gold_compile3_merged.csv. Any other table will not work!!

# OUTPUT: 1. table containing matched JGI entries for each subtype.

host_url = open(sys.argv[1],"r")
host_table = csv.reader(host_url)

host_ids = {}

for row in host_table:
	my_id = row[0].split("::")[0]
	if my_id not in host_ids:
		host_ids[my_id] = my_id

host_url.close()

host_ids = list(host_ids.keys())

jgi_url = open(sys.argv[2],"r")
jgi_table = csv.reader(jgi_url)

jgi_dict = {}
for row in jgi_table:
	jgi_dict[row[12]] = row

jgi_url.close()
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


