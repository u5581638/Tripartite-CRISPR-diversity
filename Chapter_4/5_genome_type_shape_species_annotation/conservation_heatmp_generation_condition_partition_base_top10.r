# conservation_heatmp_generation
# generate heatmap from conservation matrix for the top 10 largest clusters from leiden partitioning
# INPUT: Matrix containing conservation scores for each Pfam+DEFLOC annotated node in each partitioned cluster
#		 i.e cas12a_heatmap_matrix_top10.csv_filtered_matrix_cas12a_base3.csv
# OUTPUT: heatmap diagram showing the conservation of each protein family in each partitioned cluster for each subtype (top10 largest clusters only).

library(pheatmap)
setwd("D:/host_phage_component_partition_reconciliation/")
my_in = "cas13b_phage_only_matrix2.csv"
my_table = read.csv(my_in,row.names=1)
the_matrix <- data.matrix(my_table)
# may need to filter some columns
i = 1
# go across the rows
print(length(colnames(the_matrix)))
# Filter so that only genes in more than one partition are included!!
while (i <= length(colnames(the_matrix))) {
	switch = 0
	all_switch = 0
	k = 1

	while (k <= length(rownames(the_matrix))) {

	if (the_matrix[k,i] > 0.3 || (switch > 1 && the_matrix[k,i] > 0.15 )) { # may modify this condition to add a minimum conservation threshold
		switch = switch + 1

	}

		k = k + 1

	}
	
	if (switch < 1) {
	the_matrix <- the_matrix[,-i]
	} else {
		i = i + 1
	}
# go down the columns
	
}
print(length(colnames(the_matrix)))
print(length(rownames(the_matrix)))
# output file names
write.csv(the_matrix,paste(my_in,"filtered_matrix_cas13b_base4.csv",sep="_"))
png(paste(my_in,"_base4.png",sep="_"),res=300, units="in",width=10, height=36)

pheatmap(t(the_matrix),fontsize=16)

dev.off()

