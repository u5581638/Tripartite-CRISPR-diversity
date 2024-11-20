#Plot 4 sets of KDEs and their expected values on the same chart.
# This should then directly output the spacer distribution plot in the same general format as given by P. Fineran.

# Plan:
# 1. Compute the kernal density estimation (KDE) from the distances.
# 2. Compute the Kolmogorov-Smirnov test (KS). Plot this test result for both the sense + and sense - strand.
# 3. Colour the area undder the KDE estimations on both the positive and negative strand. Will need to manually specify the coordinates to be coloured.
# 4. Do this for complete matches as well as kmer=10,12,14.
# 4. Run this script as a batch job using R script.

library(dplyr)
library(ggplot2)

setwd("F:/spacer_expansion/non-problem_corrected_subtypes/")

graph_parameters <- function(distance_table) {
filtered_df <- filter(distance_table, mapped_strand == 1)
# need to later make negative the KDE
filtered_df2 <- filter(distance_table, mapped_strand == -1)
pos_kde <- density(filtered_df$distance)
neg_kde <- density(filtered_df2$distance)



# need to do the KS test

spacers <- distance_table %>% select(Length,distance,mapped_strand) %>% mutate(strand=factor(mapped_strand))

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




#deviation_one <- scaled_spacers %>% filter(mapped_strand==1) %>% pull(deviation)
#deviation_minus_one <- scaled_spacers %>% filter(mapped_strand==-1) %>% pull(deviation)

# The problem with this test is that it's not calculated using the original distances. Hence the difference in mapping densities isn't meaningful because the 
# sample size is arbitary.
# nearly need to compute external what the significance is, then add this to the plot.

# need to modify this by adding scaled KS result lines

# This needs to work for negative values as well. Would be good to do a KS-test for individual quadrants 
spacers <- spacers %>% filter(distance < 5000) %>% filter(distance > -5000)

n_size = length(spacers$distance)
n_size = paste("n:" , n_size, sep=" ")
#scaled_spacers <- spacers %>% mutate(distance = abs(distance)/Length) %>% select(distance,mapped_strand)

strand_one <- scaled_spacers %>% filter(mapped_strand==1) %>% pull(deviation)
strand_minus_one <- scaled_spacers %>% filter(mapped_strand==-1) %>% pull(deviation)


# scaled_spacers %>% mutate(quantile = 1-(1-distance)^2)
quantiles_strand_one <- 1-(1-strand_one)^2
quantiles_strand_minus_one <- 1-(1-strand_minus_one)^2
strand_one_ks_goodness <- ks.test(quantiles_strand_one,"punif")
strand_minus_one_ks_goodness <- ks.test(quantiles_strand_minus_one,"punif")

#strand_one_ks_goodness <- ks.test(x=pos_kde_frame$pos_kde.x,y=expected_kde_frame$expected_kde.x)
#strand_minus_one_ks_goodness <- ks.test(x=neg_kde_frame$neg_kde.x,y=neg_kde_frame$neg_kde.y)

if (strand_one_ks_goodness$p.value < 0.001) {
	sig <- "***"
} else if (strand_one_ks_goodness$p.value < 0.01) {
	sig <- "**"
} else if (strand_one_ks_goodness$p.value < 0.05) {
	sig <- "*"
} else {
	sig <- " "
}

if (strand_minus_one_ks_goodness$p.value < 0.001) {
	msig <- "***"
} else if (strand_minus_one_ks_goodness$p.value < 0.01) {
	msig <- "**"
} else if (strand_minus_one_ks_goodness$p.value < 0.05) {
	msig <- "*"
} else {
	msig <- " "
}

ks_strand_one_pre <- sig
ks_strand_minus_one_pre <- msig

# tt_strand_two_tails <- t.test((scaled_spacers %>% pull(distance))-(1/3))
upper_layer <- max(pos_kde_frame$pos_kde.y) * 0.7
bottom_layer <- min(-1 *neg_kde_frame$neg_kde.y) * 0.7
ret_out <- data.frame(
	pos_kde.x = pos_kde$x, 
	pos_kde.y = pos_kde$y,
	neg_kde.x = neg_kde$x,
	neg_kde.y = neg_kde$y,
	expected_kde.x = expected_kde$x, 
	expected_kde.y = expected_kde$y, 
	expected_minus_one_kde.x = expected_minus_one_kde$x, 
	expected_minus_one_kde.y = expected_minus_one_kde$y,
	# may need to multiply by the number of entries to form a square dataframe for the below variables
	upper_layer = upper_layer,
	bottom_layer = bottom_layer,
	ks_strand_one_pre = ks_strand_one_pre,
	ks_strand_minus_one_pre = ks_strand_minus_one_pre,
	n_size = n_size
	)
return (ret_out)
}


# complete spacer typeIIIB_cmr2.fasta
base_sequence = "typeIIIB_cmr2.fasta"

whole_sequence = paste(base_sequence, "_all_hits.csv_genomes.fasta_crisprs.lst_full_real_arr_positions.csv_distances_annotated.csv_filt_all.csv_h.csv",sep="")
whole_distance_table <- read.csv(whole_sequence)
whole_matches <- graph_parameters(whole_distance_table)
whole_match_sig_point_up.x <- which(whole_matches$pos_kde.x == max(whole_matches$pos_kde.x), arr.ind = TRUE)
whole_match_sig_point_up.y <- whole_matches$pos_kde.y[whole_match_sig_point_up.x]
whole_match_sig_point_up.x <- max(whole_matches$pos_kde.x)
whole_match_up_max <- max(whole_matches$pos_kde.y)

whole_match_down_max <- max(whole_matches$neg_kde.y)
whole_match_sig_point_down.x <- which(whole_matches$neg_kde.x == max(whole_matches$neg_kde.x), arr.ind = TRUE)
whole_match_sig_point_down.y <- whole_matches$neg_kde.y[whole_match_sig_point_down.x]
whole_match_sig_point_down.x <- max(whole_matches$neg_kde.x)

# 10 kmer matches

kmer_10_sequence = paste(base_sequence, "_detection_parallel_gadi_kmer10_31-3-2024.csv_non-self.csv_expanded.csv_corrected.csv_2_or_more_hits.csv_distances_annotated.csv_af.csv_h.csv",sep="")
kmer_10_distance_table <- read.csv(kmer_10_sequence)
kmer_10_sequence_matches <- graph_parameters(kmer_10_distance_table)
kmer_10_match_sig_point_up.x <- which(kmer_10_sequence_matches$pos_kde.x == max(kmer_10_sequence_matches$pos_kde.x), arr.ind = TRUE)
kmer_10_match_sig_point_up.y <- kmer_10_sequence_matches$pos_kde.y[kmer_10_match_sig_point_up.x]
kmer_10_match_sig_point_up.x <- max(kmer_10_sequence_matches$pos_kde.x)
kmer_10_up_max <- max(kmer_10_sequence_matches$pos_kde.y)

kmer_10_down_max <- max(kmer_10_sequence_matches$neg_kde.y)

kmer_10_match_sig_point_down.x <- which(kmer_10_sequence_matches$neg_kde.x == max(kmer_10_sequence_matches$neg_kde.x), arr.ind = TRUE)
kmer_10_match_sig_point_down.y <- kmer_10_sequence_matches$neg_kde.y[kmer_10_match_sig_point_down.x]
kmer_10_match_sig_point_down.x <- max(kmer_10_sequence_matches$neg_kde.x)

# 12 kmer matches

kmer_12_sequence = paste(base_sequence, "_detection_parallel_gadi_kmer12_31-3-2024.csv_non-self.csv_expanded.csv_corrected.csv_2_or_more_hits.csv_distances_annotated.csv_af.csv_h.csv",sep="")
kmer_12_distance_table <- read.csv(kmer_12_sequence)
kmer_12_sequence_matches <- graph_parameters(kmer_12_distance_table)
kmer_12_match_sig_point_up.x <- which(kmer_12_sequence_matches$pos_kde.x == max(kmer_12_sequence_matches$pos_kde.x), arr.ind = TRUE)
kmer_12_match_sig_point_up.y <- kmer_12_sequence_matches$pos_kde.y[kmer_12_match_sig_point_up.x]
kmer_12_match_sig_point_up.x <- max(kmer_12_sequence_matches$pos_kde.x)
kmer_12_up_max <- max(kmer_12_sequence_matches$pos_kde.y)

kmer_12_down_max <- max(kmer_12_sequence_matches$neg_kde.y)
kmer_12_match_sig_point_down.x <- which(kmer_12_sequence_matches$neg_kde.x == max(kmer_12_sequence_matches$neg_kde.x), arr.ind = TRUE)
kmer_12_match_sig_point_down.y <- kmer_12_sequence_matches$neg_kde.y[kmer_12_match_sig_point_down.x]
kmer_12_match_sig_point_down.x <- max(kmer_12_sequence_matches$neg_kde.x)

# 14 kmer matches

kmer_14_sequence = paste(base_sequence, "_detection_parallel_gadi_kmer14_31-3-2024.csv_non-self.csv_expanded.csv_corrected.csv_2_or_more_hits.csv_distances_annotated.csv_af.csv_h.csv",sep="")
kmer_14_distance_table <- read.csv(kmer_14_sequence)
kmer_14_sequence_matches <- graph_parameters(kmer_14_distance_table)
kmer_14_match_sig_point_up.x <- which(kmer_14_sequence_matches$pos_kde.x == max(kmer_14_sequence_matches$pos_kde.x), arr.ind = TRUE)
kmer_14_match_sig_point_up.y <- kmer_14_sequence_matches$pos_kde.y[kmer_14_match_sig_point_up.x]
kmer_14_match_sig_point_up.x <- max(kmer_14_sequence_matches$pos_kde.x)
kmer_14_up_max <- max(kmer_14_sequence_matches$pos_kde.y)

kmer_14_down_max <- max(kmer_14_sequence_matches$neg_kde.y)
kmer_14_match_sig_point_down.x <- which(kmer_14_sequence_matches$neg_kde.x == max(kmer_14_sequence_matches$neg_kde.x), arr.ind = TRUE)
kmer_14_match_sig_point_down.y <- kmer_14_sequence_matches$neg_kde.y[kmer_14_match_sig_point_down.x]
kmer_14_match_sig_point_down.x <- max(kmer_14_sequence_matches$neg_kde.x)

overall_max <- c(whole_match_up_max, kmer_10_up_max, kmer_12_up_max,kmer_14_up_max)
overall_max <- max(overall_max)
smallest_peak <- min(overall_max)
overall_min <- c(whole_match_down_max, kmer_10_down_max, kmer_12_down_max,kmer_14_down_max)
overall_min <- max(overall_min)
smallest_peak_neg <- min(overall_min)
whole_incre <- overall_max * 0.14 
kmer_10_incre <- whole_incre * 2
kmer_12_incre <- whole_incre * 3
kmer_14_incre <- whole_incre * 4

whole_incre_min <- overall_min * 0.14 
kmer_10_incre_neg <- whole_incre_min * 2
kmer_12_incre_neg <- whole_incre_min * 3
kmer_14_incre_neg <- whole_incre_min * 4

 


colors <- c("whole_match" = "Black","kmer=10" = "deepskyblue3","kmer=12" = "blue","kmer=14" = "violet")



deviation_plot <- ggplot(whole_matches,aes(x=pos_kde.x,y=pos_kde.y)) + theme(panel.background = element_rect(colour = "grey"),panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white"),text=element_text(size=28),axis.text.x = element_text(size = 28),axis.text.y = element_text(size = 28), legend.position="right") +
 labs(x="Distance from protospacer (bp)",y=("Protospacer density"), color="Key:") + 
 scale_color_manual(values = colors) +
#  xlim(-5000,5000) + 
  geom_vline(xintercept=0) + geom_hline(yintercept=0) +

 # geom_line(data=kmer_10_sequence_matches,aes(x=pos_kde.x,y=pos_kde.y),color="deepskyblue3",linewidth=1) + 
 # geom_area(data=kmer_10_sequence_matches,aes(x=pos_kde.x,y=pos_kde.y),fill="deepskyblue3") + 
 # geom_line(data=kmer_10_sequence_matches, aes(x=neg_kde.x,y=-neg_kde.y),color="deepskyblue3",linewidth=1) + 
#  geom_area(data=kmer_10_sequence_matches, aes(x=neg_kde.x,y=-neg_kde.y),fill="deepskyblue3") + 
#  geom_line(data=kmer_10_sequence_matches,aes(x=expected_kde.x,y=expected_kde.y, group = 1),color="deepskyblue3",linetype="dashed",linewidth=1) + 
#  geom_line(data=kmer_10_sequence_matches,aes(x=expected_minus_one_kde.x,y=-expected_minus_one_kde.y),color="deepskyblue3",linetype="dashed",linewidth=1) +
#  geom_text( aes(x=5500, y=kmer_10_match_sig_point_up.y, label=kmer_10_sequence_matches$ks_strand_one_pre[1]),size=12,color="deepskyblue3") + 
#  geom_text(aes(x=5500,y=-kmer_10_match_sig_point_down.y, label=kmer_10_sequence_matches$ks_strand_minus_one_pre[1]),size=12,color="deepskyblue3") + 
#  geom_text( aes(x=700, y=max(kmer_10_sequence_matches$pos_kde.y)*0.9, label=kmer_10_sequence_matches$n_size[1]),size=12,color="deepskyblue3") + 
#  geom_line(data=kmer_12_sequence_matches,aes(x=pos_kde.x,y=pos_kde.y),color="blue",linewidth=1) + 
#  geom_area(data=kmer_12_sequence_matches,aes(x=pos_kde.x,y=pos_kde.y),fill="blue") + 
#  geom_line(data=kmer_12_sequence_matches, aes(x=neg_kde.x,y=-neg_kde.y),color="blue",linewidth=1) + 
#  geom_area(data=kmer_12_sequence_matches, aes(x=neg_kde.x,y=-neg_kde.y),fill="blue") + 
#  geom_line(data=kmer_12_sequence_matches,aes(x=expected_kde.x,y=expected_kde.y, group = 1),color="blue",linetype="dashed",linewidth=1) + 
#  geom_line(data=kmer_12_sequence_matches,aes(x=expected_minus_one_kde.x,y=-expected_minus_one_kde.y),color="blue",linetype="dashed",linewidth=1) +
#  geom_text( aes(x=5500, y=kmer_12_match_sig_point_up.y, label=kmer_12_sequence_matches$ks_strand_one_pre[1]),size=12,color="blue") + 
#  geom_text(aes(x=5500,y=-kmer_12_match_sig_point_down.y, label=kmer_12_sequence_matches$ks_strand_minus_one_pre[1]),size=12,color="blue") + 
 # geom_text( aes(x=700, y=max(kmer_12_sequence_matches$pos_kde.y)*0.9, label=kmer_12_sequence_matches$n_size[1]),size=12,color="blue") + 


  geom_line(data=kmer_14_sequence_matches,aes(x=pos_kde.x,y=pos_kde.y),color="red",linewidth=1) + 
  geom_area(data=kmer_14_sequence_matches,aes(x=pos_kde.x,y=pos_kde.y),fill="red") + 
  geom_line(data=kmer_14_sequence_matches, aes(x=neg_kde.x,y=-neg_kde.y),color="blue2",linewidth=1) + 
  geom_area(data=kmer_14_sequence_matches, aes(x=neg_kde.x,y=-neg_kde.y),fill="blue2") + 
  geom_line(data=kmer_14_sequence_matches,aes(x=expected_kde.x,y=expected_kde.y, group = 1),color="darkred",linetype="dashed",linewidth=1) + 
  geom_line(data=kmer_14_sequence_matches,aes(x=expected_minus_one_kde.x,y=-expected_minus_one_kde.y),color="blue4",linetype="dashed",linewidth=1) +
  geom_text( aes(x=5500, y=kmer_14_match_sig_point_up.y, label=kmer_14_sequence_matches$ks_strand_one_pre[1]),size=12,color="red") + 
  geom_text(aes(x=5500,y=-kmer_14_match_sig_point_down.y, label=kmer_14_sequence_matches$ks_strand_minus_one_pre[1]),size=12,color="blue2") + 
 # geom_text( aes(x=700, y=max(kmer_14_sequence_matches$pos_kde.y)*0.9, label=kmer_14_sequence_matches$n_size[1]),size=12,color="blue4") + 



  theme(legend.position = "top")

#quantiles_strand_one <- 1-(1-strand_one)^2
#quantiles_strand_minus_one <- 1-(1-strand_minus_one)^2
#print("Coordinates:")
#print(whole_matches$upper_layer)
#print(whole_matches$upper_layer[1])
#print(whole_matches$bottom_layer[1])

# Given that this works on scaled spacers. I bet that the graphs shown in the paper actually represent scaled spacers (distance/contig length). In any case this would work in showing what we want to show!

# need to compare the expected vs abs density y values
my_out <- paste(whole_sequence,"spacer_distribution_plotdd5000_corrected_kmer14_only.png",sep="_")
ggsave(my_out,width=10,height=9,units="in")

print(kmer_10_sequence_matches$n_size)
print(kmer_12_sequence_matches$n_size)
print(kmer_14_sequence_matches$n_size)
print(whole_matches$n_size)

#strand_one_ks <- ks.test(quantiles_strand_one,"punif")
#strand_minus_one_ks <- ks.test(quantiles_strand_minus_one,"punif")


# need to plot the scaled spacer coordinates as the input to density.
# need to plot above the quantiles against the spacer coordinate density.
# Otherwise, just compute the KS score and be done with it!!

# p <- ggplot() + geom_density(data=filtered_df,aes(x=distance,y= after_stat(density), fill="red")) + geom_density(data=filtered_df2, aes(x=distance,y= after_stat(-density),fill="Cyan"))# this is only the basic version


# There are two things I need to test for statistically.
# 1. Need to test whether the distance score deviates from background based on the random probablity of two points on the same contig.
#	 a) plot observed vs simulated (do I need to simulate?) differences based on the length of the contigs.
#	 b) do KS test/T-test?
# 2. Need to test for differences between the peaks of the distribution in each quadrant.
#	a) scale the observed distances by the deviation from the expected distances. This will be obs_dist * (expected_dist-abs(observed_distance)) / expected_distance 
#	b) first compute the KDE. 
#	c) scale the points so that the maximum density is set to 1. All other points should be scaled in proportion to the maximum density. (is this really nessessary? - distorts the result!!)
#	c) create 4 seperate datasets for each quadrants.
#	d) compute the p-values based on a t-test. <- This will only be accurate if the influence of distance on the distribution is small.
 

# 3. Is there a way to incorporate both pieces of information <- correct for distance?

# Note: Protospacer density is normalised to 1. 
