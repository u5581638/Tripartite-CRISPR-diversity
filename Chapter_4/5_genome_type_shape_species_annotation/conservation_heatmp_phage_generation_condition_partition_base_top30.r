#conservation_heatmp_generation
library(pheatmap)
setwd("C:/Users/u5581/Documents/PhD/Thesis/result_chapter_2_figure_3/phage_only_virsort_heatmaps_9_9_2024/")
my_in = "casIB_new_heatmap_matrix_phage.csv"
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
write.csv(the_matrix,paste(my_in,"filtered_matrix_casIB_base2.csv",sep="_"))
png(paste(my_in,"_base.png",sep="_"),res=300, units="in",width=10, height=32)

pheatmap((the_matrix),fontsize=14)

dev.off()

