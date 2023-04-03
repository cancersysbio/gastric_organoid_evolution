#############################################################################################################################
##                                                                                                                      
##  GROWTH TRAJECTORY PLOT, ECB CULTURES
##  FOR FIGURE 5D
##                                                                                                                      
##  Date: 27 DECEMBER 2021                                                                                                                   
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

### SET PATH
setwd(file.path("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_5/Growth_trajectory_ECB/Output"))

###### FOR D2C1

### READ IN DATA
df = read.table("growth_deriv_heatmap_D2C1_Early_long_with_freq.txt",header=TRUE, sep = "\t" )

### EXTRACT READ GROUP OF INTEREST
rgs = c("RG0_R1","RG0_R2","RG0_R3","RG1_R1","RG1_R2","RG1_R3",
        "RG2_R1","RG2_R2","RG2_R3","RG3_R1","RG3_R2","RG3_R3",
        "RG4_R1","RG4_R2","RG4_R3","RG5_R1","RG5_R2","RG5_R3",
        "RG6_R1","RG6_R2","RG6_R3","RG8_R1","RG8_R2","RG8_R3") 

### SUBSET ON TIME POINTS TO INCLUDE
keep <- c("T1","T2","T3","T7","T8","T12","T13","T15","T17","T18") ### D2C1
df <- subset(df, time %in% keep)

### SET TIME AND READ GROUP REPLICATE AS FACTOR FOR PLOTTING
df$time <- factor(df$time,levels=unique(df$time))
df$RG_rep <- factor(df$RG_rep,levels=rgs)

### PLOT FIGURE
ggplot(df, aes(x=time, y=RG_rep, color=Deriv)) + 
  geom_point(aes(size=Freq)) +
  scale_colour_gradient2(low = "#0000FF", mid = "gray", high = "#FF0000",midpoint = 0.02,na.value="black") +
  scale_y_discrete(limits = rev(levels(df$RG_rep))) +

  theme(strip.text = element_text(face="bold", size=14, colour = "black",),
        strip.background = element_rect(fill="white", colour="black", size=1), 
        strip.text.y = element_text(angle = 0),
        #axis.text = element_text(colour = "black", size = 12, face = "bold" ),
        #axis.text.x = element_text(colour = "black", size = 12, face = "bold" ),
        #axis.text.y = element_text(colour = "black", size = 12, face = "bold" ),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 12, face = "bold",),
        legend.text = element_text(colour="black", size=10, face="bold"),
        panel.background = element_rect(fill = 'white'),
        #panel.grid.major = element_line(colour = "black"),
        #panel.grid.major.x = element_blank(),
        axis.line = element_line(colour = "black"),
        #axis.text.x = element_text(size=12,face = "bold"),
        legend.position="bottom") 

### SAVE FIGURE

ggsave(filename = "Dotplot_growth_trajectory_ECB_D2C1_6_5.5.pdf", width = 6, height = 5.5, units = "in", device='pdf')


###### FOR D2C2

### READ IN DATA
df = read.table("growth_deriv_heatmap_D2C2_Early_long_with_freq.txt",header=TRUE, sep = "\t" )

### EXTRACT READ GROUP OF INTEREST
rgs = c("RG0_R1","RG0_R2","RG0_R3","RG1_R1","RG1_R2","RG1_R3",
        "RG39_R1","RG39_R2","RG39_R3","RG13_R1","RG13_R2","RG13_R3",
        "RG56_R1","RG56_R2","RG56_R3","RG11_R1","RG11_R2","RG11_R3",
        "RG6_R1","RG6_R2","RG6_R3","RG9_R1","RG9_R2","RG9_R3")  ### D2C2

### SET TIME AND READ GROUP REPLICATE AS FACTOR FOR PLOTTING
df$time <- factor(df$time,levels=unique(df$time))
df$RG_rep <- factor(df$RG_rep,levels=rgs)

### PLOT FIGURE
ggplot(df, aes(x=time, y=RG_rep, color=Deriv)) + 
  geom_point(aes(size=Freq)) +
  scale_colour_gradient2(low = "#0000FF", mid = "gray", high = "#FF0000",midpoint = 0.02,na.value="black") +
  scale_y_discrete(limits = rev(levels(df$RG_rep))) +
  
  theme(strip.text = element_text(face="bold", size=14, colour = "black",),
        strip.background = element_rect(fill="white", colour="black", size=1), 
        strip.text.y = element_text(angle = 0),
        #axis.text = element_text(colour = "black", size = 12, face = "bold" ),
        #axis.text.x = element_text(colour = "black", size = 12, face = "bold" ),
        #axis.text.y = element_text(colour = "black", size = 12, face = "bold" ),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 12, face = "bold",),
        legend.text = element_text(colour="black", size=10, face="bold"),
        panel.background = element_rect(fill = 'white'),
        #panel.grid.major = element_line(colour = "black"),
        #panel.grid.major.x = element_blank(),
        axis.line = element_line(colour = "black"),
        #axis.text.x = element_text(size=12,face = "bold"),
        legend.position="bottom") 

### SAVE FIGURE

ggsave(filename = "Dotplot_growth_trajectory_ECB_D2C2_6_5.5.pdf", width = 6, height = 5.5, units = "in", device='pdf')



###### FOR D3C2

### READ IN DATA
df = read.table("growth_deriv_heatmap_D3C2_Early_long_with_freq.txt",header=TRUE, sep = "\t" )

### EXTRACT READ GROUP OF INTEREST
rgs = c("RG0_R1","RG0_R2","RG0_R3","RG1_R1","RG1_R2","RG1_R3",
        "RG2_R1","RG2_R2","RG2_R3","RG3_R1","RG3_R2","RG3_R3",
        "RG4_R1","RG4_R2","RG4_R3","RG7_R1","RG7_R2","RG7_R3",
        "RG8_R1","RG8_R2","RG8_R3")  ### D3C2

### SUBSET ON TIME POINTS TO INCLUDE
keep <- c("T1","T2","T3","T5","T6","T8","T9","T10","T13","T15","T19") ### D3C2
df <- subset(df, time %in% keep)

### SET TIME AND READ GROUP REPLICATE AS FACTOR FOR PLOTTING
df$time <- factor(df$time,levels=unique(df$time))
df$RG_rep <- factor(df$RG_rep,levels=rgs)

### PLOT FIGURE
ggplot(df, aes(x=time, y=RG_rep, color=Deriv)) + 
  geom_point(aes(size=Freq)) +
  scale_colour_gradient2(low = "#0000FF", mid = "gray", high = "#FF0000",midpoint = 0.02,na.value="black") +
  scale_y_discrete(limits = rev(levels(df$RG_rep))) +
  
  theme(strip.text = element_text(face="bold", size=14, colour = "black",),
        strip.background = element_rect(fill="white", colour="black", size=1), 
        strip.text.y = element_text(angle = 0),
        #axis.text = element_text(colour = "black", size = 12, face = "bold" ),
        #axis.text.x = element_text(colour = "black", size = 12, face = "bold" ),
        #axis.text.y = element_text(colour = "black", size = 12, face = "bold" ),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 12, face = "bold",),
        legend.text = element_text(colour="black", size=10, face="bold"),
        panel.background = element_rect(fill = 'white'),
        #panel.grid.major = element_line(colour = "black"),
        #panel.grid.major.x = element_blank(),
        axis.line = element_line(colour = "black"),
        #axis.text.x = element_text(size=12,face = "bold"),
        legend.position="bottom") 

### SAVE FIGURE

ggsave(filename = "Dotplot_growth_trajectory_ECB_D3C2_6_5.5.pdf", width = 6, height = 5.5, units = "in", device='pdf')

