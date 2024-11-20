library(pheatmap)
setwd("E:/host_phage_component_partition_reconciliation/host_partition10_final")
my_in = "cas13b_heatmap_matrix_top10_2.csv_filtered_matrix_cas13b_top10.csv"
my_table = read.csv(my_in,row.names=1)
the_matrix <- data.matrix(my_table)
png(paste(my_in,"_base10_final.png",sep="_"),res=300, units="in",width=7, height=36)

pheatmap(t(the_matrix),fontsize=17,fontsize_col=20)

dev.off()