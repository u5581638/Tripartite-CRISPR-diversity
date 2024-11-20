# function to plot expected, randomised and negative control abundances

library(dplyr)
library(ggplot2)

setwd("C:/Users/u5581/Documents/PhD/Thesis/")

base_sequence = "hex"

df <- data.frame(subtype=rep(c("Type\nV-A","Type\nV-B","Type\nV-F1","Type\nVI-B","Type\nI-B","Type\nI-D"),3),
	abundance=c(15445,1166,12642,4164,9875,3901,3146.7,107.47,3228.59,786.69,1241.43,1153.18,428,15,259,148,336,179),
	run=c(rep("Mapping",6),rep("Randomised",6),rep("Simulated",6))
	)
print(df)
ggplot(df, aes(x=subtype,y=abundance,fill=run)) + geom_bar(position="dodge", stat="identity") + theme(text = element_text(size = 20)) + geom_segment(aes(x = 0.7, y = 10000, xend = 0.7, yend = 10500)) + geom_segment(aes(x = 0.7, y = 10500, xend = 1, yend = 10500)) + geom_segment(aes(x = 1, y = 1350, xend = 1, yend = 10500)) + geom_segment(aes(x = 1.7, y = 4000, xend = 1.7, yend = 4500)) + geom_segment(aes(x = 1.7, y = 4500, xend = 2, yend = 4500)) + geom_segment(aes(x = 2, y = 1350, xend = 2, yend = 4500)) + geom_segment(aes(x = 2.7, y = 15600, xend = 2.7, yend = 16100)) + geom_segment(aes(x = 2.7, y = 16100, xend = 3, yend = 16100)) + geom_segment(aes(x = 3, y = 3200, xend = 3, yend = 16100)) + geom_segment(aes(x = 3.7, y = 1300, xend = 3.7, yend = 1800)) + geom_segment(aes(x = 3.7, y = 1800, xend = 4, yend = 1800)) + geom_segment(aes(x = 4, y = 150, xend = 4, yend = 1800)) + geom_segment(aes(x = 4.7, y = 12700, xend = 4.7, yend = 13200)) + geom_segment(aes(x = 4.7, y = 13200, xend = 5, yend = 13200)) + geom_segment(aes(x = 5, y = 3300, xend = 5, yend = 13200)) + geom_segment(aes(x = 5.7, y = 4200, xend = 5.7, yend = 4700)) + geom_segment(aes(x = 5.7, y = 4700, xend = 6, yend = 4700)) + geom_segment(aes(x = 6, y = 950, xend = 6, yend = 4700)) + geom_text( aes(x=0.8, y= 10500 + 100, label="**"),size=12) + geom_text( aes(x=1.8, y= 4500 + 100, label="**"),size=12) + geom_text( aes(x=2.8, y= 16100 + 100, label="**"),size=12) + geom_text( aes(x=3.8, y= 1800 + 100, label="**"),size=12) + geom_text( aes(x=4.8, y= 13200 + 100, label="**"),size=12) + geom_text( aes(x=5.8, y= 4700 + 100, label="**"),size=12)

my_out <- paste(base_sequence,"_custom_plotb.png",sep="_")
ggsave(my_out,width=7,height=7,units="in")