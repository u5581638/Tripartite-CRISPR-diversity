# program to go through each sorted cluster and rename cluster numbers as well as filtering out clusters which lack a single sequence greather than 300 amino acids

from Bio import SeqIO
import sys

# Helper function to rename identifiers of sequence from each cluster
def re_add (identifer):
	
	i=1
	new_id = str(identifer[0])
	while (i < len(identifer)):
		new_id = new_id + "|" + str(identifer[i])
		i += 1
	return new_id	

sequences = SeqIO.parse(sys.argv[1], "fasta")

ret_seq = []
cluster_renumber = 0 
current_cluster = "mine!"
append_switch = False 
for sequ in sequences:
	this_sequ = sequ
	this_sequ = this_sequ.description.split("|")
	clust_id = this_sequ[-1]
	this_sequ = this_sequ[:-2]
	if (current_cluster != clust_id): # new cluster
		if (len(sequ) > 300):
			append_switch = True
			cluster_renumber += 1
			print(cluster_renumber)
		else:
			print ("false")
			append_switch = False	


	if (append_switch == True):
		append_sequ = sequ
		append_sequ.description = re_add(this_sequ) + "|" + "cluster_" + str(cluster_renumber)
		ret_seq.append(append_sequ)
	current_cluster = clust_id

print("Number of sequences is: " + str(len(ret_seq)))		
SeqIO.write(ret_seq, "one_member_gt_300_renumbered_" + sys.argv[1], "fasta")


