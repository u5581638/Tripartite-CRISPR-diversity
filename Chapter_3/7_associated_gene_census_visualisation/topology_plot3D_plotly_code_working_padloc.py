import csv
import sys
import matplotlib.pyplot as plt
import plotly.express as px
import pandas as pd 

# code to generate a 3D scatter plot showing co-occurrance, abundance and distance, with each point corresponding to each representative protein sequence assigned by colour to each function category assigned by DEFLOC database searches

d3_table = pd.read_csv("0_10kb_co_occurrance_table_classified_summary_3D_input_padloc_labelled.csv")
figure = px.scatter_3d(d3_table,x='abundance',y='distance (bp)',z='co_occurrance',color='Key:',size='font_size',log_x=True,log_y=False,log_z=True)
figure.show()

#plt.savefig('3d_scatter.png', dpi=600)
#plt.show()

# setwd("C:/Users/u5581/Documents/PhD/Thesis/figure_2/figure_generation/")
# data <- read.csv("")

# mycolors <- c('red','orange','green','royalblue1','darkcyan','oldlace','purple','violet','grey','indigo')
# data$color <- mycolors[ as.numeric(data$Classification_number)]
# plot3d(x=data$abundance,y=data$co_occurrance,z=data$`distance (bp)`,
#	col=data$color,
#	type = 'p',
#	radius=.1,
#	xlab = "Abundance",ylab = "Co-occurrance",zlab="Distance")
# ggplot(co_ab,aes(x=abundance,y=co_occurrance,z=distance..bp.,color=Classification_number
# )) + theme_void() + axes_3D() + stat_3D()

# The above did not work for reasons unknown