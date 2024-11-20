
library(igraph)
library(randomcoloR)
library(dplyr)
#library(CINNA)

setwd("/Users/gaetanburgio/Documents/cluster_networks_alex_13-05-2024/cas13b/")
#edgelist <- read.csv("casIIIB_ph_filtered_16_4_2024.csv",row.names=NULL)
mygraph <- read.graph("c1.graphml",format="graphml")
# the line below 
#V(mygraph)$type <- V(mygraph)$name %in% edgelist[,1]
# detach("package:sna",unload=TRUE)
input_str = "cas13b_phage_only_partition.csv" 
out_comp <- read.csv(input_str) # read in clusters
out_comp <- out_comp %>% group_by(cluster) %>% mutate(color = randomColor()) %>% ungroup()
i = 1
while (i < length(out_comp$cluster)) {
    if (is.na(out_comp$cluster[i])) {
        out_comp$color[i] <- "#000000"
    }

    i = i + 1
}
# SET Default color if cluster == NA
# Try to increase the size of nodes
# this is the key hurdle -> need to expand the graph spacing!!
#my_weights =rep(0.8,vcount(mygraph)*2) 
mygraph <- set.edge.attribute(mygraph, "weight",index=E(mygraph),value=0.001)
my_layout <- layout_with_fr(mygraph)
my_layout <- norm_coords(my_layout, ymin=-1, ymax=1, xmin=-1, xmax=1) # normalise the coordinates to a plane between -1,1

# need to create a dataframe with the cluster numbers as column headers and 
cluster_vec <- out_comp$cluster
comp_x <- my_layout[,1]
comp_y <- my_layout[,2]
cluster_table <- data.frame(cluster_vec,comp_x,comp_y)
mean_table <- as.data.frame(aggregate(cluster_table[,2:3], list(cluster_table$cluster_vec),mean))
colnames(mean_table) = c("cluster","comp_x","comp_y")
color_table <- out_comp %>% select(cluster,color)
color_table <- unique(color_table)
mean_colors = mean_table
mean_colors <- mean_colors %>% full_join(color_table)

# need to filter here!!


ostr = paste (input_str,"color_table.csv",sep="_")
write.csv(mean_colors,ostr)
print(length(color_table))

istr = paste(input_str, "graph_phage_top_all.png",sep="_")
png(filename = istr,
    width = 2000, height = 2000, units = "px", pointsize = 5,
    bg = "white", res = NA, family = "",
    type = c("cairo"), 
    symbolfamily="default")

plot(mygraph,vertex.label.cex=0.000000002,vertex.color = out_comp$color,vertex.size=0.7,vertex.size2=0.6,edge.width=0.6,edge.arrow.size=0.07,edge.arrow.width=0.07,layout=my_layout*1,rescale=FALSE)
text(x=mean_colors$comp_x,y=mean_colors$comp_y,cex=3,labels=mean_colors$cluster,col=mean_colors$color)
legend(x=0.5,y=0.9,legend=mean_colors$cluster,col=mean_colors$color,cex=6,pch=16,ncol=3)
dev.off()
