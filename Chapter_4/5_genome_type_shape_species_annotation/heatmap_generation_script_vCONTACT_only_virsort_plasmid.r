# basic heatmap generation script
# generate heatmap from input matrix (this script is specific to the parameters used for heatmap visualisation of genome type/shape conservation for each partitioned cluster from network generation)

library(pheatmap)
setwd("D:/host_phage_batch_run_7_5_2024-27_8_2024")
my_table = read.csv("cas12a_species6_inclu_heatmap_matrix.csv",row.names=1)
#row.names(my_table) <- my_table$cluster_no
#my_table <- subset(my_table,select = -c(Component_no)) # this removes this column?
the_matrix <- data.matrix(my_table)
png(paste("cas12a_new_heatmap_matrix_host_speco6.csv","png",sep="."),res=300, units="in",width=12, height=24)
pheatmap(the_matrix,fontsize=20,cluster_rows = FALSE, cluster_cols=FALSE)

dev.off()