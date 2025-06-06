library(dplyr)
library(ggplot2)
#theme_set(theme_minimal())
#library(hrbrthemes)
#args<-commandArgs(TRUE)


# script to generate a spacer distribution. 
# Note: the parameters of the KDEs were altered for different ranges over which the distribution was measured
# INPUT: deduplicated table of PPS-spacer distances
# i.e. cas12a.fasta_whole_h.csv
# Output: Kernal density plot showing a density distribution of PPS-spacer distances centered on the PPS (the origin).

setwd("E:/spacer_expansion/final_overlay_eric_h_files_29_4_2024/")
i_sequence = "typeIA.fasta_whole_h.csv"
distance_table <- read.csv(i_sequence)
filtered_df <- filter(distance_table, mapped_strand == 1)
# need to later make negative the KDE
filtered_df2 <- filter(distance_table, mapped_strand == -1)
pos_kde <- density(filtered_df$distance) # need to use this function as the direct input
neg_kde <- density(filtered_df2$distance)

spacers <- distance_table %>% select(Length,distance,mapped_strand,Mapped_start_site,Mapped_end_site) %>% mutate(strand=factor(mapped_strand))

scaled_spacers <- spacers %>% mutate(expected_value = (distance / abs(distance)) * (Length / 3),deviation = (abs(distance) / (Length)), similarity = distance * (1 - deviation)) # %>% select(distance,mapped_strand)

pos_kde = density(spacers %>% filter(mapped_strand==1) %>% pull(distance),adjust = 1, kernel = c("gaussian"), window = kernel, n=512,from=-5000,to=5000,bw=50)
pos_kde_frame <- data.frame(pos_kde$x,pos_kde$y)
neg_kde = density(spacers %>% filter(mapped_strand==-1) %>% pull(distance), adjust = 1, kernel = c("gaussian"), window = kernel, n=512,from=-5000,to=5000,bw=50)
neg_kde_frame <- data.frame(neg_kde$x,neg_kde$y)

expected_one <- scaled_spacers %>% filter(mapped_strand==1) %>% pull(expected_value)
expected_minus_one <- scaled_spacers %>% filter(mapped_strand==-1) %>% pull(expected_value)

expected_kde <- density(expected_one, adjust = 1, kernel = c("gaussian"), window = kernel, n=512,from=-5000,to=5000, bw=50)
expected_kde_frame <- data.frame(expected_kde$x,expected_kde$y)
expected_minus_one_kde <- density(expected_minus_one,adjust = 1, kernel = c("gaussian"), window = kernel, n=512,from=-5000,to=5000, bw=50) 
expected_kde_minus_one_frame <- data.frame(expected_minus_one_kde$x,expected_minus_one_kde$y)

#filter PPS-spacer distances pairs between [-5000, 5000]????
spacers <- spacers %>% filter(distance < 5000) %>% filter(distance > -5000)
n_size = length(spacers$distance)
n_size = paste("n:" , n_size, sep=" ")

strand_one <- scaled_spacers %>% filter(mapped_strand==1) %>% pull(distance)
strand_minus_one <- scaled_spacers %>% filter(mapped_strand==-1) %>% pull(distance)
# basic mirror plot

expected_one <- scaled_spacers %>% filter(mapped_strand==1) %>% pull(expected_value)
expected_minus_one <- scaled_spacers %>% filter(mapped_strand==-1) %>% pull(expected_value)

expected_kde <- density(expected_one, adjust = 1, kernel = c("gaussian"), window = kernel, n=512,from=-5000,to=5000, bw=50)
expected_kde_frame <- data.frame(expected_kde$x,expected_kde$y)
expected_minus_one_kde <- density(expected_minus_one,adjust = 1, kernel = c("gaussian"), window = kernel, n=512,from=-5000,to=5000, bw=50) 
expected_kde_minus_one_frame <- data.frame(expected_minus_one_kde$x,expected_minus_one_kde$y)

strand_one <- scaled_spacers %>% filter(mapped_strand==1)  %>% pull(deviation)
strand_minus_one <- scaled_spacers %>% filter(mapped_strand==-1)  %>%  pull(deviation)
quantiles_strand_one <- 1-(1-strand_one)^2
quantiles_strand_minus_one <- 1-(1-strand_minus_one)^2
strand_one_ks_goodness <- ks.test(quantiles_strand_one,"punif")
strand_minus_one_ks_goodness <- ks.test(quantiles_strand_minus_one,"punif")

if (strand_one_ks_goodness$p.value < 0.001) {
	sig <- "***"
} else if (strand_one_ks_goodness$p.value < 0.01) {
	sig <- "**"
} else if (strand_one_ks_goodness$p.value < 0.05) {
	sig <- "*"
} else {
	sig <- "ns"
}

if (strand_minus_one_ks_goodness$p.value < 0.001) {
	msig <- "***"
} else if (strand_minus_one_ks_goodness$p.value < 0.01) {
	msig <- "**"
} else if (strand_minus_one_ks_goodness$p.value < 0.05) {
	msig <- "*"
} else {
	msig <- "ns"
}

ks_strand_one_pre <- paste("KS:",sig)
ks_strand_minus_one_pre <-  paste("KS:",msig)

upper_layer <- max(pos_kde_frame$pos_kde.y) * 0.7
bottom_layer <- min(-1 *neg_kde_frame$neg_kde.y) * 0.7

deviation_plot <- ggplot() + theme(panel.background = element_rect(colour = "black"),panel.grid.major=element_line(colour="grey"),panel.grid.minor=element_line(colour="grey"),text=element_text(size=28),axis.text.x = element_text(size = 28),axis.text.y = element_text(size = 28)) + labs(x="Distance from protospacer (bp)",y=("Protospacer density")) + xlim(-5000,5000) + geom_line(data=pos_kde_frame,aes(x=pos_kde.x,y=pos_kde.y),color="Black") + geom_area(data=pos_kde_frame,aes(x=pos_kde.x,y=pos_kde.y),fill="Red") + geom_line(data=neg_kde_frame, aes(x=neg_kde.x,y=-neg_kde.y),color="Black") + geom_area(data=neg_kde_frame, aes(x=neg_kde.x,y=-neg_kde.y),fill="Cyan") + geom_line(data=expected_kde_frame,aes(x=expected_kde.x,y=expected_kde.y, group = 1),color="Black") + geom_line(data=expected_kde_minus_one_frame,aes(x=expected_minus_one_kde.x,y=-expected_minus_one_kde.y),color="Black") +  geom_text( aes(x=3500, y=upper_layer, label=ks_strand_one_pre),size=12) + geom_text(aes(x=3500,y=bottom_layer, label=ks_strand_minus_one_pre),size=12) + geom_text( aes(x=3500, y=max(pos_kde_frame$pos_kde.y)*0.9, label=n_size),size=12) + theme(legend.position = "none")

quantiles_strand_one <- 1-(1-strand_one)^2
quantiles_strand_minus_one <- 1-(1-strand_minus_one)^2

# Given that this works on scaled spacers. I bet that the graphs shown in the paper actually represent scaled spacers (distance/contig length). In any case this would work in showing what we want to show!

# need to compare the expected vs abs density y values
my_out <- paste(i_sequence,"spacer_distribution_plotd_whole_5000_8_9_2024_a.png",sep="_")
ggsave(my_out,width=10,height=9,units="in")