# conservation_heatmp_generation
# script to generate a heatmap from a matrix of conservation scores
library(pheatmap)
setwd("C:/Users/u5581/Documents/PhD/Thesis/result_chapter_1_figure_3/")
my_in = "conservation_matrix2_pfam_500_all_merged.csv"
# matrix of conservation scores
my_table = read.csv(my_in,row.names=1)
the_matrix <- data.matrix(my_table)
# may need to filter some columns
i = 2
# go across the rows
print(length(colnames(the_matrix)))
while (i < length(colnames(the_matrix)) - 1) {
	switch = 0
	all_switch = 0
	k = 2
	while (k < length(rownames(the_matrix)) - 1) {
	if (the_matrix[k,i] >= 1) {
		all_switch = all_switch + 1
	}
	if (the_matrix[k,i] != 0 && the_matrix[k,i] > 0.2 && the_matrix[k,i] <= 1) { # may modify this condition to add a minimum conservation threshold
		switch = switch + 1

	}

		k = k + 1

	}
	
	if (switch < 1 || all_switch > 0) {
	the_matrix <- the_matrix[,-i]
	} else {
		i = i + 1
	}
# go down the columns
	
}
print(length(colnames(the_matrix)))
print(length(rownames(the_matrix)))
# output file names
write.csv(the_matrix,paste(my_in,"filtered_matrix_h4_twenty_percent.csv",sep="_"))
png(paste(my_in,"h4_twenty_percent.png",sep="_"),res=300, units="in",width=16, height=32)

pheatmap(t(the_matrix),fontsize=20)

dev.off()

