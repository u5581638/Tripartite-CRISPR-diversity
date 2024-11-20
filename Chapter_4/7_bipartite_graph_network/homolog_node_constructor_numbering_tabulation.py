# homolog_node_constructor
# Group spacer/phage_ids into vcontact2 identified clusters!!

import sys
import csv


with open(sys.argv[1], "r") as csvfile:
	hitmap_table = list(csv.reader(csvfile))
# .clusters host file
with open(sys.argv[2], "r") as csvfile2:
	host_net_table = list(csv.reader(csvfile2))
#.clusters phage file
with open(sys.argv[3], "r") as csvfile3:
	phage_net_table = list(csv.reader(csvfile3))

out_url = open(sys.argv[4], "w")
spam_writer = csv.writer(out_url)
spam_writer.writerow(["Genome_id","Phage_id","Genome_cluster","Phage_cluster","Phage_contig"])

out_url2 = open(sys.argv[4] + "_cluster_only.csv","w")
spam_writer2 = csv.writer(out_url2)
spam_writer2.writerow(["Genome_cluster","Phage_cluster","Phage_contig"])

# create a dictionary of host-cluster pairs. When an id is in more than one cluster, use the cluster with the lowest e-value
host_dict = {}
for host in host_net_table[1:]:
	host_cluster = "h" + host[0] # cluster number
	host_p_value = float(host[6])
	host_genomes = host[7].split(" ")
	for genome in host_genomes:
		if genome not in host_dict:
			host_dict[genome] = [host_cluster,host_p_value]
		else:
			# need to explain these p-values!!
			if (host_p_value < host_dict[genome][1]):
				host_dict[genome] = [host_cluster, host_p_value]


phage_dict = {}
for phage in phage_net_table[1:]:
	phage_cluster = "p" + phage[0]
	phage_p_value = float(phage[6])
	phage_genomes = phage[7].split(" ")
	for genome in phage_genomes:
		genome_id = genome.split(":")[0] # need to modify the table to include the phage contig_distance
		genome_id = genome_id.split("|") [1]
		if genome_id not in phage_dict:
			phage_dict[genome_id] = [phage_cluster,phage_p_value,genome] # conserving the original phage_genome_identifer.
		else:
			if (phage_p_value < phage_dict[genome_id][1]):
				phage_dict[genome_id] = [phage_cluster, phage_p_value, genome]	


# Have two cluster dictionaries. Now to lookup and write the key information to file:

for hitmap in hitmap_table[1:]:
	hitmap[0] = hitmap[0].split("|")[0]
	hitmap[1] = hitmap[1].split("|")[1]
	ret_row = [hitmap[0], hitmap[1]]
#	print("Hi")
	# This does not match to the host dict - bug!!
#	print(host_dict.keys())
#	print(hitmap[0])
	if (hitmap[0] in host_dict):
		ret_row.append(host_dict[hitmap[0]][0])
	else:
	#	print("Error: host element not found!!")
		ret_row.append(hitmap[0])
	if (hitmap[1] in phage_dict):
		ret_row.append(phage_dict[hitmap[1]][0])
		ret_row.append(phage_dict[hitmap[1]][2])
	else:
	#	print("Error: phage element not found!!")

		ret_row.append(hitmap[1])
		ret_row.append("NA")

	spam_writer.writerow(ret_row)
	spam_writer2.writerow(ret_row[2:])

out_url.close()
out_url2.close()