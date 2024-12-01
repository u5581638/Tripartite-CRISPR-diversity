library(ggplot2)
# code to generate a scatterplot comparing the effects of co-occurrance, distance and/or abundance.
co_ab <- read.csv("distance_co_occurrance.csv")
ggplot(co_ab,aes(x=distance..bp.,y=co_occurrance))+ geom_point(size=1,shape=16) + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10') + geom_density2d(bins=25) + labs(x="distance (bp)",y="co-occurrance")