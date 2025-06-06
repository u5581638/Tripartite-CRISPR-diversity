library(igraph)
#library(CINNA)

# generate a number of graph partition using the leiden clustering algorithm
# INPUT: network file in .graphml format (converted by cytoscape)
# i.e. c1.graphml
# OUTPUT: Table labelling each node in the graph by it's assigned cluster (community) as determined by leiden clustering of the network.
# i.e. cas12a_host_only_partition2.csv"
# Note: The order of entries in the table is important as this order maps to the order of nodes in the graph file. This is important later when rendering the graph coloured according to community (from leiden clustering)
setwd("D:/host_phage_component_reconciliation")
mygraph <- read.graph("c1.graphml",format="graphml")
new_graph <- as.undirected(mygraph,mode="mutual")
graph_communities <- cluster_leiden(new_graph,resolution_parameter=0.0000000001)
node <- c()
partition <- c()
i = 1
while (i <= length(graph_communities$membership)) {
	label <- graph_communities$name[i]
	val <- as.double(graph_communities$membership[i])
	node <- append(node,label)
	partition <- append(partition,val)
	i = i + 1
}
out_comp <- data.frame(node,partition)
write.csv(out_comp,"cas12a_host_only_partition2.csv")
# graph partitions should be integrated into graph as per component coloring. If possible, partitions from the same component should use the same base color.