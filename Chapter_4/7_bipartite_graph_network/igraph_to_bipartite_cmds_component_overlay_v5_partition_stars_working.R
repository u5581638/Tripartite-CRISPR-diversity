
library(igraph)
library(CINNA)
library(png)
library(magick)
library(r2r)
library(dplyr)

# script to generate and render type VI host-phage bipartite networks.
# INPUT: host-phage interaction table with cluster partition number and associated RGB color added.
# OUTPUT: rendering of host-phage interaction network using the Fruchterman Reingold layout coloured and labelled by host-phage partitioned clusters.

setwd("E:/host-phage-bipartite-coloring_11_9_2024/cas13d")
# need to create script to swap the first and second rows
input_str = "cas13b_host_phage_interaction_table.csv_names_kept2.csv_hp_colors_table_top10_only.csv"
edgelist <- read.csv(input_str,row.names=NULL)
edgelist <- select(edgelist,Phage_cluster,Genome_cluster,Phage_contig,host_cluster,host_color,phage_cluster,phage_color)
i = 1
while (i < length(edgelist$phage_cluster)) {
	if (is.na(edgelist$phage_cluster[i])) {
		edgelist$phage_color[i] <- "#000000"
	}
	i = i + 1
}
# could use just one combined table
host_colors <- edgelist$host_color
phage_colors <- edgelist$phage_color
combined_colors <- c(host_colors,phage_colors)

i = 1 
m <- hashmap()
while (i <= length(edgelist$Genome_cluster)) {
	if (!has_key(m,edgelist$Genome_cluster[i])) {
		m[[ edgelist$Genome_cluster[i] ]] <- edgelist$host_color[i]
	}
	if (!has_key(m,edgelist$Phage_cluster[i])) {
		m[[ edgelist$Phage_cluster[i] ]] <- edgelist$phage_color[i]
	}
	i = i + 1
}

host_color_list <- data.frame(edgelist$host_cluster,edgelist$host_color)
names(host_color_list) <- c("host_cluster","host_color")
host_color_list <- host_color_list %>% filter(!is.na(host_cluster))

host_color_list <- host_color_list[!duplicated(host_color_list$host_cluster), ]

phage_color_list <- data.frame(edgelist$phage_cluster,edgelist$phage_color)
names(phage_color_list) <- c("phage_cluster","phage_color")
phage_color_list <- phage_color_list %>% filter(!is.na(phage_cluster))

phage_color_list <- phage_color_list[!duplicated(phage_color_list$phage_cluster),]


new_colors <- c()
mygraph <- graph.data.frame(edgelist)
vertices <- V(mygraph)$name
i = 1 
while (i <= length(vertices)) {
	if (has_key(m,vertices[i])) {
		new_colors <- c(new_colors, m[[ vertices[i] ]])

#		print("Hello")
	}
	else {
		print("phage_color")
		print(vertices[i])
		new_colors <- c(new_colors,"#000000")
	}
	i = i + 1
}


options(max.print = .Machine$integer.max,scipen = 999)
V(mygraph)$type <- V(mygraph)$name %in% edgelist[,1]
V(mygraph)$color <- new_colors
# detach("package:sna",unload=TRUE)
# Try to increase the size of nodes
mygraph <- set.edge.attribute(mygraph, "weight",index=E(mygraph),value=0.5)
my_layout <- layout_with_fr(mygraph)
my_layout <- norm_coords(my_layout, ymin=-1, ymax=1, xmin=-1, xmax=1) # normalise the coordinates to a plane between -1,1

node_names <- V(mygraph)$name
V(mygraph)$shape <- c("rectangle","circle")[V(mygraph)$type + 1]
# url_vector <- c("cas12a.fasta_whole_h.csv_spacer_distribution_plotdd5000U.png","spcas9.fasta_whole_h.csv_spacer_distribution_plotdd5000U.png","cas12b.fasta_whole_h.csv_spacer_distribution_plotdd5000U.png","cas12f1.fasta_whole_h.csv_spacer_distribution_plotdd5000U.png","cas13b.fasta_whole_h.csv_spacer_distribution_plotdd5000U.png")

# need to create a dataframe with the component numbers as column headers and 
#my_shape_host <-shapes(shape="rect",parameters)

png(filename = paste(input_str, "host_phage_top10.png",sep="_"),
    width = 3000, height = 2000, units = "px", pointsize = 5,
    bg = "white", res = NA, family = "", restoreConsole = TRUE,
    type = c("windows", "cairo", "cairo-png"), 
    symbolfamily="default")
# replace color with 
plot(mygraph,axes=TRUE,cex=50,vertex.label.cex=0.02,vertex.color = V(mygraph)$color, vertex.size=2,vertex.size2=2 ,edge.curved = 0,edge.width=0.6,edge.arrow.size=0.07,edge.arrow.width=0.07,layout=my_layout*1,rescale=FALSE)
# host cluster legend
legend(x=0.5,y=1,legend=host_color_list$host_cluster,col=host_color_list$host_color,pch=15,cex=6,ncol=3,title="host partition no.") # check with code on the bioinfo computer
#phage cluster legend
legend(x=0.5,y=0.5,legend=phage_color_list$phage_cluster,col=phage_color_list$phage_color,pch=16,cex=6,ncol=3,title="phage partition no.")
#square circle legend()
#legend(x=-1,y=0.5,legend=c("Host","Phage"),col=c(0,1),pch=c(15,16))
# need to save in 3000*3000 plot dimensions during export to see the plot properly. Still needs to be larger (experiment with xlim/ylim) but is a signficant improvement!

#text(x=mean_table[,2],y=mean_table[,3],cex=10,labels=mean_table[,1])
dev.off()


