# basic heatmap generation script
library(pheatmap)
setwd("D:/cas12b_heatmap_rerun_15_4_2024/heatmap_repeat/phage")
my_table = read.csv("out_components.csv",row.names=NULL)
#row.names(my_table) <- my_table$cluster_no
my_table <- subset(my_table,select = -c(Component_no)) # this removes this column?
the_matrix <- data.matrix(my_table)
png("out_component_matrix.png",res=300, units="in",width=16, height=32)
pheatmap(the_matrix)

dev.off()