# script to extract a single representative protein sequence from a given cluster
from Bio import SeqIO
import sys
import csv
import copy
import statistics
import math

# function to extract the most representative sequence in each cluster, based on the smallest eucidilean distances derived from the pairwise similarities between sequences in each cluster.
def rep_sequence_determinant(mds_table_url, maximum_cluster_bound=0.06, maximum_cluster_bound_x=None, maximum_cluster_bound_y=None):
	with open(mds_table_url, "r") as csvfile:
		hit_table = list(csv.reader(csvfile)) # input is a csv file with sequence_id, MDS x coordinate and MDS y coordinate on each row.

	hit_table = hit_table[1:]
	point_dict = {}
	for point in hit_table:
		hit_table_cpy = copy.deepcopy(hit_table)
		point_id = point[0]
		point_x_coord = point[1]
		point_y_coord = point[2]
		value_list = []
		for point_cpy in hit_table_cpy:
			if (point_cpy[0] != point_id): # this line may not end up being nessessary
				x_distance = float(point_x_coord) - float(point_cpy[1])
				y_distance =  float(point_y_coord) - float(point_cpy[2])
				euc_dis = math.sqrt((x_distance)**2 + (y_distance)**2)
				value_list.append([point_id, point_cpy[0], euc_dis, x_distance, y_distance, point_cpy[1], point_cpy[2], point_x_coord, point_y_coord, point_id]) # may want to add sequence x, y coords??
			else:
				print("self matching hit!!")
		point_dict[point_id] = value_list # the value list should be a list of lists. Essentially a 2D matrix for each sequence.
	print("Checkpoint!!")
	# now assemble points into clusters		
	ret_list = []
	rm_list = set()
	node_num = 0
	other_num = 0
	for point in point_dict: 
		if (point in rm_list):
			continue
		nodes = {point} 
		edges = point_dict[point] 
		stack = []
		stack = copy.deepcopy(edges)
		switch = 0
		print("start main loop")
		while(len(stack) != 0):
			aux_node = stack.pop() # pop the last element on the stack
			if (aux_node[1] not in nodes and aux_node[1] not in rm_list and float(aux_node[2]) < maximum_cluster_bound): 
				nodes.add(aux_node[1])
				node_num += 1
			
				if (aux_node[1] in point_dict and aux_node[1] not in rm_list):
					stack.extend(point_dict[aux_node[1]])
					switch = 1
				else:
					print("missing node, removed previously??")
			else:
				other_num += 1

		ret_list.append(list(nodes))
		print("progress!")
		for ele in nodes:
			if (ele not in rm_list):
				rm_list.add(ele) 
	print("pre_sort!!")
	ret_list.sort(key=len, reverse=True)
	# substitute in original MDS table entries!!
	print("Sorted!!")
	i = 0
	while (i < len(ret_list)):
		k = 0
		while (k < len(ret_list[i])):
			vari = point_dict[ret_list[i][k]]
			vari = copy.deepcopy(vari)
			ret_list[i][k] =  vari
			k += 1
		i += 1
	print("pre_clustering!!")
	ret_cluster_avg_ids = []
	for cluster in ret_list:
		avg_x_y_ele = []
		for table in cluster:
			x_ele = float(table[0][-3])
			y_ele = float(table[0][-2])
			id_ele = table[0][-1]
			avg_x_y_ele.append([id_ele,x_ele,y_ele])
		avg_x_coord = statistics.mean(list(map(lambda x: float(x[1]), avg_x_y_ele)))
		avg_y_coord = statistics.mean(list(map(lambda x: float(x[1]), avg_x_y_ele)))
		distance_euc = {} 
		for ele in avg_x_y_ele:
			distance_euc[math.sqrt((avg_x_coord - ele[1])**2 + (avg_y_coord - ele[2])**2)] = id_ele
		lowest_distance = min(distance_euc.keys())
		lowest_distance_id = distance_euc[lowest_distance]
		ret_cluster_avg_ids.append(lowest_distance_id)
	return ret_cluster_avg_ids	
			




		







