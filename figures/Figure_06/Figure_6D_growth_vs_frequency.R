#############################################################################################################################
##                                                                                                                      
##  COMPARISON SUBCLONE GROWTH AND FREQUENCY
##  FOR FIGURE 6D, ED9D, ED10C, ED11C
##                                                                                                                      
##  Date: 27 DECEMBER 2021                                                                                                                   
##  
##  Author: Kasper Karlsson
##
##                                                                                                                      
############################################################################################################################


library(stats)
library(reshape2)
library(ggplot2)
library(forcats)
library(heatmap3)
library("gplots")
library(ggpubr)
library(gridExtra)
library(grid)
library(scales)


setwd("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_6/Figure_6D/norm_RGC_matrices")


### READ IN COLORS (COMMON ALL PLOTS)

colors <- read.csv("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/sequtils/Supporting_files/barcode_colors_hexcode_RG_match.csv",header=FALSE)
colnames(colors) <- c("RG","color")
colors$RG <- paste0("RG",colors$RG)


############## FIG 6D - D2C2R2 - GROWTH vs. FREQ TOP GROWING CLONES

sample <- "D2C2"
freq_data_input <- read.table("D2C2_RGC_matrix_norm.txt",header=TRUE,as.is=TRUE) 
freq_data_filter <- freq_data_input[which(freq_data_input$D2C2_Parent1>1e-5),] 
freq_data_filter$RG <- rownames(freq_data_filter)

### ADD CORRECT COLORS
freq_data_merge <- merge(freq_data_filter,colors,by="RG")
freq_data <- freq_data_merge[order(match(freq_data_merge$RG,freq_data_filter$RG)),]

### TIME POINTS USED
parent <- paste0(sample,"_Parent1")
R1T1 <- paste0(sample,"_R1_T1")
R2T1 <- paste0(sample,"_R2_T1")
R3T1 <- paste0(sample,"_R3_T1")
R1T2 <- paste0(sample,"_R1_T2")
R2T2 <- paste0(sample,"_R2_T2")
R3T2 <- paste0(sample,"_R3_T2")
R1T4 <- paste0(sample,"_R1_T4")
R2T4 <- paste0(sample,"_R2_T4")
R3T4 <- paste0(sample,"_R3_T4")
R1T5 <- paste0(sample,"_R1_T5")
R2T5 <- paste0(sample,"_R2_T5")
R3T5 <- paste0(sample,"_R3_T5")
R1T8 <- paste0(sample,"_R1_T8")
R2T8 <- paste0(sample,"_R2_T8")
R3T8 <- paste0(sample,"_R3_T8")
R1T10 <- paste0(sample,"_R1_T10")
R2T10 <- paste0(sample,"_R2_T10")
R3T10 <- paste0(sample,"_R3_T10")
R1T12 <- paste0(sample,"_R1_T12")
R2T12 <- paste0(sample,"_R2_T12")
R3T12 <- paste0(sample,"_R3_T12")
R1T19 <- paste0(sample,"_R1_T19")
R2T19 <- paste0(sample,"_R2_T19")
R3T19 <- paste0(sample,"_R3_T19")
R1T20 <- paste0(sample,"_R1_T20")
R2T20 <- paste0(sample,"_R2_T20")
R3T20 <- paste0(sample,"_R3_T20")

### FUNTION TO PLOT FIGURE

plot_gr_freq_rg_color <- function(freq_data,freq_tp,previous_tp,name){
  freq_data_loop <- freq_data
  freq_data_loop["growth"] <- freq_data[,freq_tp] / freq_data[,previous_tp] # CALCULATE GROWTH (BETWEEN THE TWO SPECIFIED TIME POINTS)
  tg <- c("RG0","RG1","RG6","RG9","RG11","RG13","RG39","RG56") # SPECIFY HIGHLIGHTED SUBCLONES
  topGrowers <-subset(freq_data_loop, RG %in% tg) # SUBSET HIGHLIGHTED SUBCLONES
  topGrowers_RG0 <- subset(freq_data_loop, RG=="RG0") # CAN BE USED TO CHANGE PARAMETER FOR A SPECIFIC RG
  ### PLOT
  p <- ggplot()
  p=p+geom_point(data=freq_data_loop, aes(x=freq_data_loop[,freq_tp], y=freq_data_loop$growth),size=0.5, color="gray") # ALL RG SMALLER SIZE AND GRAY
  p=p+geom_point(data=topGrowers,aes(x=topGrowers[,freq_tp], y=topGrowers$growth),color=topGrowers$color,size=2) # TOP PLOT RG SEPARATELY, LARGER SIZE
  p=p+geom_point(data=topGrowers_RG0,aes(x=topGrowers_RG0[,freq_tp], y=topGrowers_RG0$growth),color=topGrowers_RG0$color,size=2) # PLOT RG0 SEPARATELY, IF NEEDED
  p=p+scale_x_continuous(trans = 'log10',
                         breaks = trans_breaks('log10', function(x) 10^x),
                         labels = trans_format('log10', math_format(10^.x)))
  p=p+scale_y_continuous(trans = 'log10',
                         breaks = trans_breaks('log10', function(x) 10^x),
                         labels = trans_format('log10', math_format(10^.x)))
  
  p=p+ggtitle(name) 
  p=p+theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
            panel.background = element_rect(fill = 'white'))
  
  return(p)
}

### SAVE FIGURE
pdf(paste0("Growth_Freq_Specific_RG_D2C2R2T1_gray_3_3.pdf"),width=3,height=3,pointsize=0.1)
plot_gr_freq_rg_color(freq_data, R2T1, parent,"R2T1")
dev.off()


############## FIG ED 9D  - D2C1R1 - GROWTH vs. FREQ TOP GROWING CLONES

sample <- "D2C1"
freq_data_input <- read.table("D2C1_RGC_matrix_norm.txt",header=TRUE,as.is=TRUE) 
freq_data_filter <- freq_data_input[which(freq_data_input$D2C1_Parent1>1e-5),] 
freq_data_filter$RG <- rownames(freq_data_filter)

### ADD CORRECT COLORS
freq_data_merge <- merge(freq_data_filter,colors,by="RG")
freq_data <- freq_data_merge[order(match(freq_data_merge$RG,freq_data_filter$RG)),]

### TIME POINTS USED
parent <- paste0(sample,"_Parent1")
R1T1 <- paste0(sample,"_R1_T1")
R2T1 <- paste0(sample,"_R2_T1")
R3T1 <- paste0(sample,"_R3_T1")
R1T2 <- paste0(sample,"_R1_T2")
R2T2 <- paste0(sample,"_R2_T2")
R3T2 <- paste0(sample,"_R3_T2")
R1T4 <- paste0(sample,"_R1_T4")
R2T4 <- paste0(sample,"_R2_T4")
R3T4 <- paste0(sample,"_R3_T4")
R1T5 <- paste0(sample,"_R1_T5")
R2T5 <- paste0(sample,"_R2_T5")
R3T5 <- paste0(sample,"_R3_T5")
R1T8 <- paste0(sample,"_R1_T8")
R2T8 <- paste0(sample,"_R2_T8")
R3T8 <- paste0(sample,"_R3_T8")
R1T10 <- paste0(sample,"_R1_T10")
R2T10 <- paste0(sample,"_R2_T10")
R3T10 <- paste0(sample,"_R3_T10")
R1T12 <- paste0(sample,"_R1_T12")
R2T12 <- paste0(sample,"_R2_T12")
R3T12 <- paste0(sample,"_R3_T12")
R1T19 <- paste0(sample,"_R1_T19")
R2T19 <- paste0(sample,"_R2_T19")
R3T19 <- paste0(sample,"_R3_T19")
R1T20 <- paste0(sample,"_R1_T20")
R2T20 <- paste0(sample,"_R2_T20")
R3T20 <- paste0(sample,"_R3_T20")

### PLOT FIGURE

plot_gr_freq_rg_color <- function(freq_data,freq_tp,previous_tp,name){
  freq_data_loop <- freq_data
  freq_data_loop["growth"] <- freq_data[,freq_tp] / freq_data[,previous_tp]
  tg <- c("RG0","RG1","RG2","RG4") 
  topGrowers <-subset(freq_data_loop, RG %in% tg)
  topGrowers_RG0 <- subset(freq_data_loop, RG=="RG0") 
  print (topGrowers)
  p <- ggplot()
  p=p+geom_point(data=freq_data_loop, aes(x=freq_data_loop[,freq_tp], y=freq_data_loop$growth),size=0.5, color="gray") 
  p=p+geom_point(data=topGrowers,aes(x=topGrowers[,freq_tp], y=topGrowers$growth),color=topGrowers$color,size=2)
  p=p+geom_point(data=topGrowers_RG0,aes(x=topGrowers_RG0[,freq_tp], y=topGrowers_RG0$growth),color=topGrowers_RG0$color,size=2)
  p=p+scale_x_continuous(trans = 'log10',
                         breaks = trans_breaks('log10', function(x) 10^x),
                         labels = trans_format('log10', math_format(10^.x)))
  p=p+scale_y_continuous(trans = 'log10',
                         breaks = trans_breaks('log10', function(x) 10^x),
                         labels = trans_format('log10', math_format(10^.x)))
  
  p=p+ggtitle(name) 
  p=p+theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
            panel.background = element_rect(fill = 'white'))
  
  return(p)
}

### SAVE FIGURE
pdf(paste0("Growth_Freq_Specific_RG_D2C1R1T1_gray_3_3.pdf"),width=3,height=3,pointsize=0.1)
plot_gr_freq_rg_color(freq_data, R1T1, parent,"R1T1")
dev.off()


############## FIG ED 10C - D2C1R2 - GROWTH vs. FREQ TOP GROWING CLONES

sample <- "D2C1"
freq_data_input <- read.table("D2C1_RGC_matrix_norm.txt",header=TRUE,as.is=TRUE) 
freq_data_filter <- freq_data_input[which(freq_data_input$D2C1_Parent1>1e-5),] 
freq_data_filter$RG <- rownames(freq_data_filter)

### ADD CORRECT COLORS
freq_data_merge <- merge(freq_data_filter,colors,by="RG")
freq_data <- freq_data_merge[order(match(freq_data_merge$RG,freq_data_filter$RG)),]

### TIME POINTS USED
parent <- paste0(sample,"_Parent1")
R1T1 <- paste0(sample,"_R1_T1")
R2T1 <- paste0(sample,"_R2_T1")
R3T1 <- paste0(sample,"_R3_T1")
R1T2 <- paste0(sample,"_R1_T2")
R2T2 <- paste0(sample,"_R2_T2")
R3T2 <- paste0(sample,"_R3_T2")
R1T4 <- paste0(sample,"_R1_T4")
R2T4 <- paste0(sample,"_R2_T4")
R3T4 <- paste0(sample,"_R3_T4")
R1T5 <- paste0(sample,"_R1_T5")
R2T5 <- paste0(sample,"_R2_T5")
R3T5 <- paste0(sample,"_R3_T5")
R1T8 <- paste0(sample,"_R1_T8")
R2T8 <- paste0(sample,"_R2_T8")
R3T8 <- paste0(sample,"_R3_T8")
R1T10 <- paste0(sample,"_R1_T10")
R2T10 <- paste0(sample,"_R2_T10")
R3T10 <- paste0(sample,"_R3_T10")
R1T12 <- paste0(sample,"_R1_T12")
R2T12 <- paste0(sample,"_R2_T12")
R3T12 <- paste0(sample,"_R3_T12")
R1T19 <- paste0(sample,"_R1_T19")
R2T19 <- paste0(sample,"_R2_T19")
R3T19 <- paste0(sample,"_R3_T19")
R1T20 <- paste0(sample,"_R1_T20")
R2T20 <- paste0(sample,"_R2_T20")
R3T20 <- paste0(sample,"_R3_T20")

### PLOT FIGURE

plot_gr_freq_rg_color <- function(freq_data,freq_tp,previous_tp,name){
  freq_data_loop <- freq_data
  freq_data_loop["growth"] <- freq_data[,freq_tp] / freq_data[,previous_tp]
  tg <- c("RG0","RG1","RG2","RG4") #D2C1
  topGrowers <-subset(freq_data_loop, RG %in% tg)
  topGrowers_RG0 <- subset(freq_data_loop, RG=="RG0") 
  print (topGrowers)
  p <- ggplot()
  p=p+geom_point(data=freq_data_loop, aes(x=freq_data_loop[,freq_tp], y=freq_data_loop$growth),size=0.5, color="gray") 
  p=p+geom_point(data=topGrowers,aes(x=topGrowers[,freq_tp], y=topGrowers$growth),color=topGrowers$color,size=2)
  p=p+geom_point(data=topGrowers_RG0,aes(x=topGrowers_RG0[,freq_tp], y=topGrowers_RG0$growth),color=topGrowers_RG0$color,size=2)
  p=p+scale_x_continuous(trans = 'log10',
                         breaks = trans_breaks('log10', function(x) 10^x),
                         labels = trans_format('log10', math_format(10^.x)))
  p=p+scale_y_continuous(trans = 'log10',
                         breaks = trans_breaks('log10', function(x) 10^x),
                         labels = trans_format('log10', math_format(10^.x)))
  
  p=p+ggtitle(name) 
  p=p+theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
            panel.background = element_rect(fill = 'white'))
  
  return(p)
}

### SAVE FIGURE
pdf(paste0("Growth_Freq_Specific_RG_D2C1R2T1_gray_3_3.pdf"),width=3,height=3,pointsize=0.1)
plot_gr_freq_rg_color(freq_data, R2T1, parent,"R2T1")
dev.off()


############## FIG ED 11C - D3C2R1 - GROWTH vs. FREQ TOP GROWING CLONES

sample <- "D3C2"
freq_data_input <- read.table("D3C2_RGC_matrix_norm.txt",header=TRUE,as.is=TRUE) 
freq_data_filter <- freq_data_input[which(freq_data_input$D3C2_Parent>1e-5),]
freq_data_filter$RG <- rownames(freq_data_filter)

### ADD CORRECT COLORS
freq_data_merge <- merge(freq_data_filter,colors,by="RG")
freq_data <- freq_data_merge[order(match(freq_data_merge$RG,freq_data_filter$RG)),]


### TIME POINTS USED
parent <- paste0(sample,"_Parent")
R1T1 <- paste0(sample,"_R1_T1")
R2T1 <- paste0(sample,"_R2_T1")
R3T1 <- paste0(sample,"_R3_T1")
R1T2 <- paste0(sample,"_R1_T2")
R2T2 <- paste0(sample,"_R2_T2")
R3T2 <- paste0(sample,"_R3_T2")
R1T4 <- paste0(sample,"_R1_T4")
R2T4 <- paste0(sample,"_R2_T4")
R3T4 <- paste0(sample,"_R3_T4")
R1T5 <- paste0(sample,"_R1_T5")
R2T5 <- paste0(sample,"_R2_T5")
R3T5 <- paste0(sample,"_R3_T5")
R1T8 <- paste0(sample,"_R1_T8")
R2T8 <- paste0(sample,"_R2_T8")
R3T8 <- paste0(sample,"_R3_T8")
R1T10 <- paste0(sample,"_R1_T10")
R2T10 <- paste0(sample,"_R2_T10")
R3T10 <- paste0(sample,"_R3_T10")
R1T12 <- paste0(sample,"_R1_T12")
R2T12 <- paste0(sample,"_R2_T12")
R3T12 <- paste0(sample,"_R3_T12")
R1T19 <- paste0(sample,"_R1_T19")
R2T19 <- paste0(sample,"_R2_T19")
R3T19 <- paste0(sample,"_R3_T19")
R1T20 <- paste0(sample,"_R1_T20")
R2T20 <- paste0(sample,"_R2_T20")
R3T20 <- paste0(sample,"_R3_T20")

### PLOT FIGURE

plot_gr_freq_rg_color <- function(freq_data,freq_tp,previous_tp,name){
  freq_data_loop <- freq_data
  freq_data_loop["growth"] <- freq_data[,freq_tp] / freq_data[,previous_tp]
  tg <- c("RG0","RG1","RG2","RG3","RG4","RG7","RG8")
  topGrowers <-subset(freq_data_loop, RG %in% tg) 
  topGrowers_RG0 <- subset(freq_data_loop, RG=="RG0") 
  print (topGrowers)
  p <- ggplot()
  p=p+geom_point(data=freq_data_loop, aes(x=freq_data_loop[,freq_tp], y=freq_data_loop$growth),size=0.5, color="gray") 
  p=p+geom_point(data=topGrowers,aes(x=topGrowers[,freq_tp], y=topGrowers$growth),color=topGrowers$color,size=2)
  p=p+geom_point(data=topGrowers_RG0,aes(x=topGrowers_RG0[,freq_tp], y=topGrowers_RG0$growth),color=topGrowers_RG0$color,size=2)
  p=p+scale_x_continuous(trans = 'log10',
                         breaks = trans_breaks('log10', function(x) 10^x),
                         labels = trans_format('log10', math_format(10^.x)))
  p=p+scale_y_continuous(trans = 'log10',
                         breaks = trans_breaks('log10', function(x) 10^x),
                         labels = trans_format('log10', math_format(10^.x)))
  
  p=p+ggtitle(name) 
  p=p+theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
            panel.background = element_rect(fill = 'white'))
  
  return(p)
}

### SAVE FIGURE
pdf(paste0("Growth_Freq_Specific_RG_D3C2R1T1_gray_3_3.pdf"),width=3,height=3,pointsize=0.1)
plot_gr_freq_rg_color(freq_data, R1T1, parent,"R1T1")
dev.off()
