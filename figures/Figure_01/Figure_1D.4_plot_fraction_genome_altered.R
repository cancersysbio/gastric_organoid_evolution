#############################################################################################################################
##                                                                                                                      
##  PLOT FRACTION GENOME ALTERED
##                                                                                                                      
##  Date: 27 DECEMBER 2021                                                                                                                   
##  
##  Author: Kasper Karlsson
##
##                                                                                                                      
############################################################################################################################


library("RColorBrewer")
library(ggplot2)
library(dplyr)
library(scales)

setwd(file.path("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_1/PLOTS/Figure_1D_Fraction_genome_altered"))
df <- read.table("Summary_fraction_genome_altered_per_win_movingAve_25p.csv",header=TRUE,comment.char = "",sep=",")

### REMOVE ECB AND WT SAMPLES
df <- subset(df,ECB_sample != "YES")
df <- subset(df,Clone != "WT")

df$Day <- as.numeric(df$Day) ### MAKE DAYS NUMERIC
color = c("D1"="#00B050","D2"="#843C0C","D3"="#7C7C7C") ### ADD CLONE COLORS

### PLOT FIGURE
ggplot(df, aes(x=Day, y=fraction_total_win_altered, color=Donor, shape=Clone))+
  geom_line() +
  geom_point(size = 3) +
  scale_color_manual(values=color) +
  scale_x_continuous(limits = c(0, 900)) +
  theme(strip.text = element_text(face="bold", size=14, colour = "black",),
        strip.background = element_rect(fill="white", colour="black", size=1), 
        axis.text = element_text(colour = "black", size = 12, face = "bold" ),
        axis.title = element_text(colour = "black", size = 16, face = "bold" ),
        axis.title.y =  element_blank(),
        axis.title.x =  element_blank(),
        plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 12, face = "bold",),
        legend.text = element_text(colour="black", size=10, face="bold"),
        panel.background = element_rect(fill = 'white'),
        axis.line = element_line(colour = "black"),
        legend.position="bottom") 

### SAVE FIGURE
ggsave(filename = "Figure1D_fraction_genome_altered_6_3.pdf", width = 6, height = 3, units = "in", device='pdf')
