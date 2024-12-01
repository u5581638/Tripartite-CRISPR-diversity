# generate binomial plot and assess significance by quadrant.
 
library(dplyr)
library(ggplot2)
library(hrbrthemes)
library(ggpubr)
library(ggprism)
library(EMT)


setwd("D:/spacer_expansion/finished_problem_corrected_subtypes/finished_2")
i_sequence <- "spcas9_detection_parallel_gadi_kmer14_29-3-2024.csv_non-self.csv_expanded.csv_corrected.csv_2_or_more_hits.csv_justified.csv_distances_annotated.csv_af.csv_h.csv"
distance_table <- read.csv(i_sequence)
distance_table <- distance_table %>% filter(distance > -5000) %>% filter(distance < 5000)
spacers <- distance_table %>% select(Length,distance,mapped_strand) %>% mutate(strand=factor(mapped_strand))
quandrant_upper <- spacers %>% filter(mapped_strand == 1)
quandrant_lower <- spacers %>% filter(mapped_strand == -1)
total_size = length(spacers$distance)
quandrant_upper_size <- length(quandrant_upper$distance)
quandrant_lower_size <- length(quandrant_lower$distance)
coord_position <- max(c(quandrant_upper_size,quandrant_lower_size)) + 1

# Bionomial probability
# P(X) = ((n!) /(n-x)!x!)*pxqn-x
strand_probability <- binom.test(quandrant_lower_size,total_size,1/2) 

# need to draw and save a strand probability plot



# compute multinomial significance!!
# apply a similar procedure to label each plot!!
# may want to limit to [-50000,50000]

# quad_probability <- multinomial.test(c(quadrant_1_size,quadrant_2_size,quadrant_3_size,quadrant_4_size),c(0.25,0.25,0.25,0.25), MonteCarlo=TRUE,ntrial=500000)

# now convert to *


# now need to add the binomial calculated p-values
 
if (strand_probability$p.value < 0.001) {
     pastreik <- "***"
 } else if (strand_probability$p.value < 0.01){
     pastreik <- "**"
 } else if (strand_probability$p.value < 0.05){
     pastreik <- "*"
 } else {
    pastreik <- "ns"
}
#pastreik <- paste(signif(strand_probability$p.value,3), pastreik,sep = " ")
quandrant_strands <- data.frame(strand_bias = c("Non-target","Target"), abundance = c( quandrant_upper_size,quandrant_lower_size))
p_frame <- data.frame(Positive = c(quandrant_lower_size),Negative = c(quandrant_upper_size), p_value = c(pastreik), y.position=c(coord_position))
# need to color plot and label with bionomial significance!!

strand_bar_chart <- ggplot(quandrant_strands,aes(strand_bias,abundance)) + geom_col(fill="Black") + geom_bracket(xmin = 1,xmax = 2,y.position=coord_position *1.1 + 3,label=pastreik) + theme(text=element_text(size=20))

my_out <- paste(i_sequence,"target_distribution_14kmer_plot.png",sep="_")
ggsave(my_out,width=4,height=7,units="in")


# compute multinomial significance!!
# apply a similar procedure to label each plot!!
# may want to limit to [-50000,50000]

quadrant_1 <- spacers %>% filter(mapped_strand == 1 & distance > 0 & distance < 5000) # counter_clockwise from the positive distance- positive strand quadrant
quadrant_1_size <- length(quadrant_1$distance)
quadrant_2 <- spacers %>% filter(mapped_strand == 1 & distance < 0 & distance > -5000)
quadrant_2_size <- length(quadrant_2$distance)
quadrant_3 <- spacers %>% filter(mapped_strand == -1 & distance < 0 & distance > -5000)
quadrant_3_size <- length(quadrant_3$distance)
quadrant_4 <- spacers %>% filter(mapped_strand == -1 & distance > 0 & distance < 5000)
quadrant_4_size <- length(quadrant_4$distance)


# quad_probability <- multinomial.test(c(quadrant_1_size,quadrant_2_size,quadrant_3_size,quadrant_4_size),c(0.25,0.25,0.25,0.25), MonteCarlo=TRUE,ntrial=500000)


label_pos <- max(c(quadrant_1_size,quadrant_2_size)) * 1.2 + 30
label_neg <- -1 * max(c(quadrant_3_size, quadrant_4_size)) * 1.2 - 30
quads_pos <- data.frame(strand = c("Positive","Negative"), abundance = c( quadrant_1_size,quadrant_2_size))
quads_neg <- data.frame(strand = c("Positive","Negative"), abundance = c( quadrant_4_size,quadrant_3_size))

# need to compute bionomial probabilities for each quadrant combination
total_si <- quadrant_1_size + quadrant_2_size + quadrant_3_size + quadrant_4_size
quad_1_p <- binom.test(quadrant_1_size, total_si,1/4)
quad_2_p <- binom.test(quadrant_2_size, total_si,1/4)
quad_3_p <- binom.test(quadrant_3_size, total_si,1/4)
quad_4_p <- binom.test(quadrant_4_size, total_si,1/4)


# now convert to *

if (quad_1_p$p.value < 0.001) {
	qp1 <- "***"
} else if (quad_1_p$p.value < 0.01) {
	qp1 <- "**"
} else if (quad_1_p$p.value < 0.05) {
	qp1 <- "*"
} else {
	qp1 <- "ns"
}

if (quad_2_p$p.value < 0.001) {
	qp2 <- "***"
} else if (quad_2_p$p.value < 0.01) {
	qp2 <- "**"
} else if (quad_2_p$p.value < 0.05) {
	qp2 <- "*"
} else {
	qp2 <- "ns"
}

if (quad_3_p$p.value < 0.001) {
	qp3 <- "***"
} else if (quad_3_p$p.value < 0.01) {
	qp3 <- "**"
} else if (quad_3_p$p.value < 0.05) {
	qp3 <- "*"
} else {
	qp3 <- "ns"
}

if (quad_4_p$p.value < 0.001) {
	qp4 <- "***"
} else if (quad_4_p$p.value < 0.01) {
	qp4 <- "**"
} else if (quad_4_p$p.value < 0.05) {
	qp4 <- "*"
} else {
	qp4 <- "ns"
}

# need to draw C-bars -> has to be done manually because 

#quad_2_3_h1_df <- data.frame(x1 = 0.5,x2 = 0.7,y1 = -1 * quadrant_3_size, y2 = -1 * quadrant_3_size) 
#quad_2_3_v1_df <- data.frame(x1 = 0.5,x2 = 0.7,y1 = 1 * quadrant_2_size, y2 = 1 * quadrant_2_size) 
#quad_2_3_h1_df <- data.frame(x1=0.5,x2=0.5,y1=quadrant_2_size,y2=quadrant_3_size)

print(quadrant_1_size)
print(quadrant_2_size)
print(label_pos)
print(quadrant_3_size) 
print(quadrant_4_size)

# draw the base plot:
quad_bar_strand <- ggplot() + theme(text=element_text(size=20),axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18)) + 
geom_col(data=quads_neg, aes(x=strand,y=-abundance,fill="Red")) + 
geom_col(data=quads_pos, aes(x=strand,y=abundance, fill="Blue")) +
scale_y_continuous("no. protospacers") + geom_label( aes(x=1, y=label_pos, label='5\'')) + geom_label( aes(x=2, y=label_pos, label='3\'')) + geom_label( aes(x=1, y=label_neg, label='3\'')) + geom_label( aes(x=2, y=label_neg, label='5\'')) + geom_hline(yintercept=0) 

quad_bar_strand <- quad_bar_strand + 
geom_text( aes(x=2, y=quadrant_1_size * 1.1 + 5, label=qp1),size=10) + geom_text( aes(x=1, y=quadrant_2_size * 1.1 + 5, label=qp2),size=10) + geom_text( aes(x=2, y=-quadrant_4_size  * 1.1 - 5, label=qp4),size=10) + geom_text( aes(x=1, y=-quadrant_3_size * 1.1 - 5, label=qp3),size=10)

my_out <- paste(i_sequence,"quad_distribution_14kmer_plot5000.png",sep="_")
ggsave(my_out,width=5,height=7,units="in")

# now need to add the binomial calculated p-values
 
