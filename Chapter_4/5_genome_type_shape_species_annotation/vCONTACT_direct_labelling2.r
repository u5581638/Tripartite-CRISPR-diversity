# vCONTACT_direct_labelling
library(randomcoloR)
library(dplyr)
library(igraph)
library(CINNA)

setwd("D:/cas12b_heatmap_rerun_15_4_2024/heatmap_repeat/phage/")
host_cluster_graph = read.graph("cas12b_c1_graph.graphml",format="graphml")
host_virsort_table = read.csv("final_table.csv")
host_plasme_table = read.csv("cas12b.fasta_faidx_bp_window.fasta_plasmid_overlaps.fa_p.csv")
host_gold_annotation_table = read.csv("cas12b_phage_c1.clusters_gold_annotations.csv")
host_ncbi_annotation_table = read.csv("cas12b_phage_c1.clusters_ncbi_annotations.csv")
colnames(host_gold_annotation_table) <- c("genome_id","genome","sample_info","sra","genbank","gproject","annotation","category","gsample","gbank","unnamed1","unnamed2","unnamed3","unnamed4","unnamed5","unnamed6","unnamed7","unnamed8","unnamed9","unnamed10","unnamed11","unnamed12","unnamed13","unnamed14","unnamed15","location","sample_description","sample_environment","sample_summary","sample_source","sample_creation","sample_culture","oxygen","sample_substance","biome","metagenome_type","metagenome_substance","type","creation2","sample_area")
colnames(host_ncbi_annotation_table) <- c("genome_id","genbank_proj","bioproject_id","biosample_id","genbank_id","na","abundance1","abundance2","annotation","strain","plasmid","version","completeness","majority","fraction","date","version_id","institute","GCF_no","matches","ftp_url","frameshift","type_material","date2","polypoidy","genome_type","number1","number2","identity","zero_column","number3","number4","storage_base","Uni","date3","number5","number6","number7","number8")


main_annotation_frame = data.frame(matrix(nrow = 0, ncol = 3))
colnames(main_annotation_frame) <- c("vert","genome_id","annotation")
graph_vertexes <- V(host_cluster_graph)
graph_vertexes[1]
i = 1
while (i <= length(graph_vertexes)) {
	my_vertex <- names(graph_vertexes[i])
	print(my_vertex)
	vertex_id <- unlist(strsplit(my_vertex,split = ":"))
	vertex_id <- vertex_id[1]
	print("Hello")
	if (vertex_id %in% host_gold_annotation_table$genome_id) {
		print("if1")
		print(which(host_gold_annotation_table$genome_id == vertex_id))
		anno_row_index <- which(host_gold_annotation_table$genome_id == vertex_id)
	#	anno_row <- host_gold_annotation_table[anno_row_index]
		print("if1.2")
		anno_row_species <- host_gold_annotation_table$unnamed6[anno_row_index]
		anno_row_sample <- host_gold_annotation_table$sample_environment[anno_row_index]
		print(anno_row_species)
		print(anno_row_sample)
		if (is.na(anno_row_species) && ! is.na (anno_row_sample)) {
			anno_row <- anno_row_sample
		} else if (! is.na(anno_row_species) && is.na(anno_row_sample) ) {
			anno_row <- anno_row_species
		} else {
			anno_row = ""
		} 
		new_frame <- data.frame(vert=c(my_vertex),genome_id=c(vertex_id),annotation=c(anno_row))
		main_annotation_frame <- rbind(main_annotation_frame,new_frame)	
		} else if (vertex_id %in% host_ncbi_annotation_table$genome_id){
		print("if2")
		anno_row_index <- which(host_ncbi_annotation_table$genome_id == vertex_id)
		#anno_row_ <- host_ncbi_annotation_table[anno_row_index]
		anno_row_species <- host_ncbi_annotation_table$annotation[anno_row_index]
		anno_row_sample <- host_ncbi_annotation_table$genome_type[anno_row_index]
		if (is.na(anno_row_species) && ! is.na(anno_row_sample) ) {
			anno_row <- anno_row_sample
		} else if (! is.na(anno_row_species) && is.na(anno_row_sample) ) {
			anno_row <- anno_row_species
		} else {
			anno_row = ""
		}
		new_frame <- data.frame(vert=c(my_vertex),genome_id=c(vertex_id),annotation=c(anno_row))
		main_annotation_frame <- rbind(main_annotation_frame, new_frame)
	}
	else {
		print("if3")
		new_frame <- data.frame(vert=c(my_vertex),genome_id=c(vertex_id),annotation=c(""))
		main_annotation_frame <- rbind(main_annotation_frame, new_frame)
	}
	print("Anno_frame:")

	i = i + 1
}
#key is whether the group keeps the order synchronised. If this doesn't work the table must be resorted
main_annotation_frame <- main_annotation_frame %>% group_by(annotation) %>% mutate(color = randomColor()) %>% ungroup()

#vertex.color = main_annotation_frame$color
# plot(host_cluster_graph,vertex.label.cex=0.4,vertex.color = main_annotation_frame$color)
plot(host_cluster_graph,vertex.label.cex=0.08,,vertex.size=4,vertex.color = main_annotation_frame$color,edge.width=0.005,layout=layout.circle,main="circle",edge.arrow.size=0.005,edge.arrow.width=0.005,rescale=TRUE)
# exSinglesOnly <- 
#    ex %>% 
#    group_by(id,day) %>% # the complete group of interest
#    mutate(duplicate = n()) %>% # count number in each group
#    filter(duplicate == 1) %>% # select only unique records
#    select(-duplicate) # remove group count column
color_frame <- main_annotation_frame %>% select(annotation,color)
color_frame <- color_frame %>% distinct(annotation,color)
legend("topleft",legend=color_frame$annotation,pch=1,col=color_frame$color,cex=0.62
)
