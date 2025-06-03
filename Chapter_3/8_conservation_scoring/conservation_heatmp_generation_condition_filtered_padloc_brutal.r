library(pheatmap)
setwd("D:/host_phage_component_partition_reconciliation/")
# render heatmap of conservation matrix with tight conservation thresholds
# INPUT: conservation matrix with at leat one well greater than 20% conserved
# i.e.
my_in = "cas12a_combined_matrix.csv"
my_table = read.csv(my_in,row.names=1)
the_matrix <- data.matrix(my_table)
# OUTPUT: heatmap rendering in PNG format of strongly filtered conservation matrix.
# i.e
png(paste(my_in,"h4_filtered_brutal.png",sep="_"),res=300, units="in",width=16, height=32)
# generate a heatmap from a matrix of filtered conservation scores
pheatmap(t(the_matrix),fontsize=20)

dev.off()