# basic heatmap generation script

# generate heatmap from input matrix
# INPUT: conservation matrix
# i.e. see below
# OUTPUT: heatmap showing partitioned community specific information (generic script)
library(pheatmap)
setwd("D:/host_phage_batch_run_7_5_2024/interaction_host_graphs")
my_in = "cas12a_comp_matrix_host.csv"
my_table = read.csv(my_in,row.names=NULL)
#row.names(my_table) <- my_table$cluster_no
my_table <- subset(my_table,select = -c(cluster_no)) # this removes this column?
the_matrix <- data.matrix(my_table)
png(paste(my_in,"h.png",sep="_"),res=300, units="in",width=8, height=32)
pheatmap(the_matrix)

dev.off()