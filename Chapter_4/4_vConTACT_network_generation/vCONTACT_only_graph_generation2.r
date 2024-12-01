library(igraph)
library(randomcoloR)
library(dplyr)
#library(CINNA)

# render networks using the Fruchterman Reingold layout
setwd("E:/host_phage_batch_run_7_5_2024/host_only_vCONTACT_networks/casIIIA")
#edgelist <- read.csv("cas12b_ph_filtered_16_4_2024.csv",row.names=NULL)
mygraph <- read.graph("c1.graphml",format="graphml")
# the line below 
#V(mygraph)$type <- V(mygraph)$name %in% edgelist[,1]
# detach("package:sna",unload=TRUE)
input_str = "casIIIA_host_only_components.csv" 
out_comp <- read.csv(input_str) # read in components
out_comp <- out_comp %>% group_by(component) %>% mutate(color = randomColor()) %>% ungroup()
# Try to increase the size of nodes
# this is the key hurdle -> need to expand the graph spacing!!
#my_weights =rep(0.8,vcount(mygraph)*2) 
mygraph <- set.edge.attribute(mygraph, "weight",index=E(mygraph),value=0.001)
my_layout <- layout_with_fr(mygraph)
my_layout <- norm_coords(my_layout, ymin=-1, ymax=1, xmin=-1, xmax=1) # normalise the coordinates to a plane between -1,1

# need to create a dataframe with the component numbers as column headers and 
component_vec <- out_comp$component
comp_x <- my_layout[,1]
comp_y <- my_layout[,2]
component_table <- data.frame(component_vec,comp_x,comp_y)
mean_table <- as.data.frame(aggregate(component_table[,2:3], list(component_table$component_vec),mean))
colnames(mean_table) = c("component","comp_x","comp_y")
color_table <- out_comp %>% select(component,color)
color_table <- unique(color_table)
mean_colors = mean_table
mean_colors <- mean_colors %>% full_join(color_table)
print(length(color_table))

istr = paste(input_str, "graph_exp3.png",sep="_")
png(filename = istr,
    width = 2000, height = 2000, units = "px", pointsize = 5,
    bg = "white", res = NA, family = "",
    type = c("cairo"), 
    symbolfamily="default")

plot(mygraph,vertex.label.cex=0.000000002,vertex.color = out_comp$color,vertex.size=0.7,vertex.size2=0.6,edge.width=0.6,edge.arrow.size=0.07,edge.arrow.width=0.07,layout=my_layout*1,rescale=FALSE)
text(x=mean_colors$comp_x,y=mean_colors$comp_y,cex=10,labels=mean_colors$component,col=mean_colors$color)
dev.off()
