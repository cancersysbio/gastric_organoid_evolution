#############################################################################################################################
##                                                                                                                      
##  GROWTH TRAJECTORY PLOT, NON-ECB CULTURES
##  FOR FIGURES 3B AND S15A
##                                                                                                                      
##  
##  Author: Kasper Karlsson
##
##                                                                                                                      
############################################################################################################################


library(ggplot2)
library(reshape2)
library(heatmap3)
library(ggpubr)
library(tidyverse)
library(tidyr)


setwd(file.path("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_3/Figure_3B_Growth_trajectory_nonECB/Final_Output"))

df = read.csv("nonECB_growth_deriv_heatmap_EarlyLate_raw_counts.csv",header=TRUE) # Figure 3B
df = read.csv("nonECB_growth_deriv_heatmap_EarlyMid_raw_counts.csv",header=TRUE) # Figure S9A 


df <- df[order(df$oriPass),]
df <- df[order(df$sample),]

df2 <- subset(df, growth_FC != "Na")
df2$growth_FC <- as.numeric(df2$growth_FC)

#### ORDER
unique(df2$sample)

df2$sample <- factor(df2$sample, levels = c("D3C3","D3C2","D3C1","D2C3","D2C2","D2C1","D1C3","D1C2","D1C1"))


### PLOT

ggplot(df2, aes(x=newPass, y=sample, color=growthDeriv)) + 
  geom_point(aes(size=growth_FC))+
  scale_colour_gradient2(low = "#0000FF", mid = "gray", high = "#FF0000",midpoint = 0.075,na.value="black")+
  scale_size(limits = c(0,150), breaks = c(0,50,100,150)) +
  

  theme(strip.text = element_text(face="bold", size=14, colour = "black",),
        strip.background = element_rect(fill="white", colour="black", size=1), 
        axis.text = element_text(colour = "black", size = 12, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 12, face = "bold" ),
        axis.text.y = element_text(colour = "black", size = 12, face = "bold" ),
        axis.title = element_blank(),
        plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 12, face = "bold",),
        legend.text = element_text(colour="black", size=10, face="bold"),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(colour = "black")) 

ggsave(filename = "Dotplot_growth_trajectory_nonECB_growthFC_EarlyLate_Adjusted_8_4.pdf", width = 8, height = 4, dpi = 300, units = "in", device='pdf')
ggsave(filename = "Dotplot_growth_trajectory_nonECB_growthFC_EarlyMid_4_4.pdf", width = 4, height = 4, dpi = 300, units = "in", device='pdf')


