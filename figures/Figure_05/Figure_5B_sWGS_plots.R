#############################################################################################################################
##                                                                                                                      
##  GENERATE SCATTERPLOTS FROM SHALLOW WGS RESULTS PROCESSED WITH QDNASEQ
##                                                                                                                      
##  Date: 15 AUGUST 2021                                                                                                                   
##  
##  Author: Moritz Przybilla
##
##                                                                                                                      
############################################################################################################################

# clear workspace
rm(list=ls())
set.seed(16011985) # set the seed of random number generator 

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("reshape2", "optparse", "BSgenome", "RColorBrewer", "ggplot2", "scales", "DescTools", "dendextend", "tidyverse", 
                      "Matrix", "devtools", "Matrix.utils", "matrixStats", "readr", "magrittr", "BiocManager", 
                      "biomaRt", "httr", "ComplexHeatmap")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages, repos = "http://cran.us.r-project.org")
lapply(list.of.packages, require, character.only = TRUE)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

# add function
`%notin%` <- Negate(`%in%`)


#####################################################################################
# READ IN DATA OF INTEREST FROM SHALLOW WGS
#####################################################################################
# set input directory and collect the files
# shallow
wgs.files <- list.files("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_5/INPUT/CNA_smooth_50kb_shallow", pattern = "CNA_smooth.txt", recursive = F, full.names = T, all.files = T)
all.wgs.files <- wgs.files
o.dir <- "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_5/PLOTS/Figure_5B_CNV/"
dir.create(paste0(o.dir, "qDNAseq_plots_updated"))


#####################################################################################
# FIGURE 5B D2C1 R1 - R3, T3, T8, T20
#####################################################################################



sample.ids <- c("D2C1_R1_189d", "D2C1_R2_189d",  "D2C1_R3_189d",
                "D2C1_R1_258d", "D2C1_R2_258d",  "D2C1_R3_258d",
                "D2C1_R1_427d", "D2C1_R2_427d",  "D2C1_R3_427d")  
# add a sample id column
# read the data in
tables <- all.wgs.files[grep(paste(sample.ids[1:9], collapse = "|"), all.wgs.files)]
replicate <- str_split_fixed(basename(tables), "_", 6)[,2]
timepoint <- str_split_fixed(basename(tables), "_", 6)[,3]
tables <- lapply(tables, read.table, stringsAsFactors = F, header = T)
complete.wgs.list <- list()

replicate
timepoint
tables

j <- 1
for (j in 1:length(tables)){
  
  table <- tables[[j]]
  colnames(table)[ncol(table)] <- "cnVal"
  table$replicate <- replicate[j]
  table$timepoint <- timepoint[j]
  print(paste0(replicate[j], "_", timepoint[j]))
  
  # transform the data to get approx cnv states
  wgs.cnVal<- 2*2**table$cnVal
  wgs.cn.states <- wgs.cnVal
  
  # set up a dataframe with all information needed
  wgs.data <- data.frame(cn_state=wgs.cn.states, chr=table$chromosome, coordinates=c(1:50497), replicate = table$replicate, timepoint = table$timepoint)
  
  # Upperlimit the data to 8 copies
  wgs.data$cn_state[which(wgs.data$cn_state>8)] <- 8
  wgs.data$cn_state[which(wgs.data$cn_state < 0)] <- 0
  
  #This is needed to plot the chromosomes in the right order
  wgs.data$order <- c(1:50497)
  
  complete.wgs.list[[j]] <- wgs.data
}

complete.wgs.data <- bind_rows(complete.wgs.list)

timepoint

replicates <- c("R1", "R2", "R3")
replicate <- replicates[1]
for (replicate in replicates){
  replicate.data <- complete.wgs.data[grep(replicate, complete.wgs.data$replicate),]
  replicate.data <- replicate.data[!is.na(replicate.data$replicate), ]
  replicate.data$timepoint <- factor(replicate.data$timepoint, levels = c("189d", "258d", "427d"))
  
  # create a ggplot for this
  plotbulk<- ggplot(replicate.data, aes(coordinates, cn_state)) +
    geom_point(aes(colour = cn_state), size = 0.01) +
    scale_colour_gradient2(low = "midnightblue", high = "deeppink4", mid ="#F8F8F8", midpoint = 2, limits=c(0,4)) +
    theme_bw() +
    scale_y_continuous(breaks = c(0, 2, 4), limits = c(0,4), position = "left") +
    theme(panel.grid.major = element_blank(),
          # strip.text = element_text(face="bold", size=20, colour = "black"),
          strip.text = element_blank(),
          strip.background = element_rect(fill="white", colour="black", size=1), 
          axis.title = element_text(colour = "black", size = 16, face = "bold" ),
          plot.title = element_text(colour = "black", size = 12, face = "bold", hjust = 0.5),
          legend.title = element_text(color = "black", size = 9, face = "bold",),
          legend.text = element_text(colour="black", size=7, face="bold"),
          legend.position = "None",
          axis.text = element_text(colour = "black", size = 14, face = "bold" ),
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.spacing.x=unit(0, "lines"), panel.border = element_rect(linetype =3))+
    facet_grid(timepoint ~ reorder(chr,order), scales="free", space="free") + 
    ggExtra::removeGrid() + labs(x="Chromosomes", y="", color = "Copy Number")+
    scale_x_continuous(expand = c(0.01, 0.01))
  
  pdf(paste0(o.dir, "qDNAseq_plots_updated/WGS_D2C1_", replicate, "_WGS_Scatterplot.pdf"), width = 5.5, height = 12.5)
  print(plotbulk)
  dev.off()
  
  ### FOR SUPPLEMENTARY FIGURE 16
  plotbulk<- ggplot(replicate.data, aes(coordinates, cn_state)) +
    geom_point(aes(colour = cn_state), size = 0.01) +
    scale_colour_gradient2(low = "midnightblue", high = "deeppink4", mid ="#F8F8F8", midpoint = 2, limits=c(0,4)) +
    theme_bw() +
    scale_y_continuous(breaks = c(0, 2, 4), limits = c(0,4), position = "right") +
    theme(panel.grid.major = element_blank(),
          strip.text = element_text(face="bold", size=24, colour = "black",),
          strip.background = element_rect(fill="white", colour="black", size=1),
          strip.text.x = element_blank(),
          axis.text.y = element_text(colour = "black", size = 18, face = "bold" ),
          legend.text = element_text(colour="black", size=7, face="bold"),
          legend.position = "None",
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.spacing.x=unit(0, "lines"), panel.border = element_rect(linetype =3))+
    facet_grid(timepoint ~ reorder(chr,order), scales="free", space="free", switch="both") +
    ggExtra::removeGrid() + labs(x="", y="", color = "Copy Number")+
    scale_x_continuous(expand = c(0.01, 0.01))
  
  pdf(paste0("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_5/PLOTS/Figure_S16_CNV/qDNAseq_plots_updated/", "WGS_D2C1_", replicate,"_S16_WGS_Scatterplot_10_5.pdf"), width = 10, height = 5)
  print(plotbulk)
  dev.off()
  
}

#####################################################################################
# FIGURE 5B D2C2 R1 - R3 T3, T8, T12
#####################################################################################


sample.ids <- c("D2C2_R1_189d", "D2C2_R2_189d",  "D2C2_R3_189d",
                "D2C2_R1_258d", "D2C2_R2_258d",  "D2C2_R3_258d",
                "D2C2_R1_315d", "D2C2_R2_315d",  "D2C2_R3_315d")  
# add a sample id column
# read the data in
tables <- all.wgs.files[grep(paste(sample.ids[1:9], collapse = "|"), all.wgs.files)]
replicate <- str_split_fixed(basename(tables), "_", 6)[,2]
timepoint <- str_split_fixed(basename(tables), "_", 6)[,3]
tables <- lapply(tables, read.table, stringsAsFactors = F, header = T)
complete.wgs.list <- list()

replicate
timepoint
tables

j <- 1
for (j in 1:length(tables)){
  
  table <- tables[[j]]
  colnames(table)[ncol(table)] <- "cnVal"
  table$replicate <- replicate[j]
  table$timepoint <- timepoint[j]
  print(paste0(replicate[j], "_", timepoint[j]))
  
  # transform the data to get approx cnv states
  wgs.cnVal<- 2*2**table$cnVal
  wgs.cn.states <- wgs.cnVal
  
  # set up a dataframe with all information needed
  wgs.data <- data.frame(cn_state=wgs.cn.states, chr=table$chromosome, coordinates=c(1:50497), replicate = table$replicate, timepoint = table$timepoint)
  
  # Upperlimit the data to 8 copies
  wgs.data$cn_state[which(wgs.data$cn_state>8)] <- 8
  wgs.data$cn_state[which(wgs.data$cn_state < 0)] <- 0
  
  #This is needed to plot the chromosomes in the right order
  wgs.data$order <- c(1:50497)
  
  complete.wgs.list[[j]] <- wgs.data
}

complete.wgs.data <- bind_rows(complete.wgs.list)

replicates <- c("R1", "R2", "R3")
replicate <- replicates[1]
for (replicate in replicates){
  
  replicate.data <- complete.wgs.data[grep(replicate, complete.wgs.data$replicate),]
  replicate.data <- replicate.data[!is.na(replicate.data$replicate), ]
  replicate.data$timepoint <- factor(replicate.data$timepoint, levels = c("189d", "258d", "315d"))
  
  # create a ggplot for this
  plotbulk<- ggplot(replicate.data, aes(coordinates, cn_state)) +
    geom_point(aes(colour = cn_state), size = 0.01) +
    scale_colour_gradient2(low = "midnightblue", high = "deeppink4", mid ="#F8F8F8", midpoint = 2, limits=c(0,4)) +
    theme_bw() +
    scale_y_continuous(breaks = c(0, 2, 4), limits = c(0,4), position = "left") +
    theme(panel.grid.major = element_blank(),
          # strip.text = element_text(face="bold", size=20, colour = "black"),
          strip.text = element_blank(),
          strip.background = element_rect(fill="white", colour="black", size=1), 
          axis.title = element_text(colour = "black", size = 16, face = "bold" ),
          plot.title = element_text(colour = "black", size = 12, face = "bold", hjust = 0.5),
          legend.title = element_text(color = "black", size = 9, face = "bold",),
          legend.text = element_text(colour="black", size=7, face="bold"),
          legend.position = "None",
          axis.text = element_text(colour = "black", size = 14, face = "bold" ),
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.spacing.x=unit(0, "lines"), panel.border = element_rect(linetype =3))+
    facet_grid(timepoint ~ reorder(chr,order), scales="free", space="free") + 
    ggExtra::removeGrid() + labs(x="Chromosomes", y="", color = "Copy Number")+
    scale_x_continuous(expand = c(0.01, 0.01))
  
  pdf(paste0(o.dir, "qDNAseq_plots_updated/WGS_D2C2_", replicate, "_WGS_Scatterplot.pdf"), width = 5.5, height = 12.5)
  print(plotbulk)
  dev.off()
  
  ### FOR SUPPLEMENTARY FIGURE 16
  plotbulk<- ggplot(replicate.data, aes(coordinates, cn_state)) +
    geom_point(aes(colour = cn_state), size = 0.01) +
    scale_colour_gradient2(low = "midnightblue", high = "deeppink4", mid ="#F8F8F8", midpoint = 2, limits=c(0,4)) +
    theme_bw() +
    scale_y_continuous(breaks = c(0, 2, 4), limits = c(0,4), position = "right") +
    theme(panel.grid.major = element_blank(),
          strip.text = element_text(face="bold", size=24, colour = "black",),
          strip.background = element_rect(fill="white", colour="black", size=1),
          strip.text.x = element_blank(),
          axis.text.y = element_text(colour = "black", size = 18, face = "bold" ),
          legend.text = element_text(colour="black", size=7, face="bold"),
          legend.position = "None",
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.spacing.x=unit(0, "lines"), panel.border = element_rect(linetype =3))+
    facet_grid(timepoint ~ reorder(chr,order), scales="free", space="free", switch="both") +
    ggExtra::removeGrid() + labs(x="", y="", color = "Copy Number")+
    scale_x_continuous(expand = c(0.01, 0.01))
  
  pdf(paste0("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_5/PLOTS/Figure_S16_CNV/qDNAseq_plots_updated/", "WGS_D2C2_", replicate,"_S16_WGS_Scatterplot_10_5.pdf"), width = 10, height = 5)
  print(plotbulk)
  dev.off()
  
}



#####################################################################################
# FIGURE 5B D3C2 R1 - R3 T3, T8, T12
#####################################################################################

sample.ids <- c("D3C2_R1_189d", "D3C2_R2_189d",  "D3C2_R3_189d",
                "D3C2_R1_259d", "D3C2_R2_259d",  "D3C2_R3_259d",
                "D3C2_R1_441d", "D3C2_R2_329d",  "D3C2_R3_441d")  
# add a sample id column
# read the data in
tables <- all.wgs.files[grep(paste(sample.ids[1:9], collapse = "|"), all.wgs.files)]
replicate <- str_split_fixed(basename(tables), "_", 6)[,2]
timepoint <- str_split_fixed(basename(tables), "_", 6)[,3]
tables <- lapply(tables, read.table, stringsAsFactors = F, header = T)
complete.wgs.list <- list()

replicate
timepoint
tables

j <- 1
for (j in 1:length(tables)){
  
  table <- tables[[j]]
  colnames(table)[ncol(table)] <- "cnVal"
  table$replicate <- replicate[j]
  table$timepoint <- timepoint[j]
  print(paste0(replicate[j], "_", timepoint[j]))
  
  # transform the data to get approx cnv states
  wgs.cnVal<- 2*2**table$cnVal
  wgs.cn.states <- wgs.cnVal
  
  # set up a dataframe with all information needed
  wgs.data <- data.frame(cn_state=wgs.cn.states, chr=table$chromosome, coordinates=c(1:50497), replicate = table$replicate, timepoint = table$timepoint)
  
  # Upperlimit the data to 8 copies
  wgs.data$cn_state[which(wgs.data$cn_state>8)] <- 8
  wgs.data$cn_state[which(wgs.data$cn_state < 0)] <- 0
  
  #This is needed to plot the chromosomes in the right order
  wgs.data$order <- c(1:50497)
  
  complete.wgs.list[[j]] <- wgs.data
}

complete.wgs.data <- bind_rows(complete.wgs.list)

replicates <- c("R1", "R2", "R3")
replicate <- replicates[1]
for (replicate in replicates){
  
  replicate.data <- complete.wgs.data[grep(replicate, complete.wgs.data$replicate),]
  replicate.data <- replicate.data[!is.na(replicate.data$replicate), ]
  replicate.data$timepoint <- factor(replicate.data$timepoint, levels = c("189d", "259d", "329d", "441d"))
  
  # create a ggplot for this
  plotbulk<- ggplot(replicate.data, aes(coordinates, cn_state)) +
    geom_point(aes(colour = cn_state), size = 0.01) +
    scale_colour_gradient2(low = "midnightblue", high = "deeppink4", mid ="#F8F8F8", midpoint = 2, limits=c(0,4)) +
    theme_bw() +
    scale_y_continuous(breaks = c(0, 2, 4), limits = c(0,4), position = "left") +
    theme(panel.grid.major = element_blank(),
          # strip.text = element_text(face="bold", size=20, colour = "black"),
          strip.text = element_blank(),
          strip.background = element_rect(fill="white", colour="black", size=1), 
          axis.title = element_text(colour = "black", size = 16, face = "bold" ),
          plot.title = element_text(colour = "black", size = 12, face = "bold", hjust = 0.5),
          legend.title = element_text(color = "black", size = 9, face = "bold",),
          legend.text = element_text(colour="black", size=7, face="bold"),
          legend.position = "None",
          axis.text = element_text(colour = "black", size = 14, face = "bold" ),
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.spacing.x=unit(0, "lines"), panel.border = element_rect(linetype =3))+
    facet_grid(timepoint ~ reorder(chr,order), scales="free", space="free") + 
    ggExtra::removeGrid() + labs(x="Chromosomes", y="", color = "Copy Number")+
    scale_x_continuous(expand = c(0.01, 0.01))
  
  pdf(paste0(o.dir, "qDNAseq_plots_updated/WGS_D3C2_", replicate, "_WGS_Scatterplot.pdf"), width = 5, height = 12.5)
  print(plotbulk)
  dev.off()
  
  ### FOR SUPPLEMENTARY FIGURE 16
  plotbulk<- ggplot(replicate.data, aes(coordinates, cn_state)) +
    geom_point(aes(colour = cn_state), size = 0.01) +
    scale_colour_gradient2(low = "midnightblue", high = "deeppink4", mid ="#F8F8F8", midpoint = 2, limits=c(0,4)) +
    theme_bw() +
    scale_y_continuous(breaks = c(0, 2, 4), limits = c(0,4), position = "right") +
    theme(panel.grid.major = element_blank(),
          strip.text = element_text(face="bold", size=24, colour = "black",),
          strip.background = element_rect(fill="white", colour="black", size=1),
          strip.text.x = element_blank(),
          axis.text.y = element_text(colour = "black", size = 18, face = "bold" ),
          legend.text = element_text(colour="black", size=7, face="bold"),
          legend.position = "None",
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.spacing.x=unit(0, "lines"), panel.border = element_rect(linetype =3))+
    facet_grid(timepoint ~ reorder(chr,order), scales="free", space="free", switch="both") +
    ggExtra::removeGrid() + labs(x="", y="", color = "Copy Number")+
    scale_x_continuous(expand = c(0.01, 0.01))
  
  pdf(paste0("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_5/PLOTS/Figure_S16_CNV/qDNAseq_plots_updated/", "WGS_D3C2_", replicate,"_S16_WGS_Scatterplot_10_5.pdf"), width = 10, height = 5)
  print(plotbulk)
  dev.off()
  
}

#####################################################################################
# FIGURE 5B TOP CNV PLOTS
#####################################################################################


sample.ids <- c("D2C1_Parent","D2C2_Parent","D3C2_Parent")


i <- 1
for (i in 1:length(sample.ids)){
  
  # add a sample id column
  # read the data in
  tables <- all.wgs.files[grep(sample.ids[i], all.wgs.files)]
  sample.tmp <- sample.ids[i]
  print (tables)
  tables <- lapply(tables, read.table, stringsAsFactors = F, header = T)
  complete.wgs.list <- list()
  print (length(tables))
  j <- 1
  for (j in 1:length(tables)){
    
    table <- tables[[j]]
    colnames(table)[ncol(table)] <- "cnVal"
    table$timepoint <- sample.tmp
    #print(sample.tmp)
    
    # transform the data to get approx cnv states
    wgs.cnVal<- 2*2**table$cnVal
    wgs.cn.states <- wgs.cnVal
    
    # set up a dataframe with all information needed
    wgs.data <- data.frame(cn_state=wgs.cn.states, chr=table$chromosome, coordinates=c(1:50497), timepoint = table$timepoint)
    
    # Upperlimit the data to 8 copies
    wgs.data$cn_state[which(wgs.data$cn_state>8)] <- 8
    wgs.data$cn_state[which(wgs.data$cn_state < 0)] <- 0
    
    #This is needed to plot the chromosomes in the right order
    wgs.data$order <- c(1:50497)
    #print (wgs.data)
    complete.wgs.list[[j]] <- wgs.data
  }
  
  #print (complete.wgs.list)
  complete.wgs.data <- bind_rows(complete.wgs.list)
  
  # create a ggplot for this
  plotbulk<- ggplot(complete.wgs.data, aes(coordinates, cn_state)) +
    geom_point(aes(colour = cn_state), size = 0.01) +
    scale_colour_gradient2(low = "midnightblue", high = "deeppink4", mid ="#F8F8F8", midpoint = 2, limits=c(0,4)) +
    theme_bw() +
    scale_y_continuous(breaks = c(0, 2, 4), limits = c(0,4), position = "left") +
    theme(panel.grid.major = element_blank(),
          strip.text = element_text(face="bold", size=0, colour = "black",),
          strip.background = element_rect(fill="white", colour="white", size=1), 
          axis.text = element_text(colour = "black", size = 14, face = "bold" ),
          axis.title = element_text(colour = "black", size = 12, face = "bold" ),
          plot.title = element_text(colour = "black", size = 12, face = "bold", hjust = 0.5),
          legend.title = element_text(color = "black", size = 9, face = "bold",),
          legend.text = element_text(colour="black", size=7, face="bold"),
          legend.position = "None",
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.spacing.x=unit(0, "lines"), panel.border = element_rect(linetype =3))+
    facet_grid( ~ reorder(chr,order), scales="free", space="free", switch="both") + 
    ggExtra::removeGrid() + labs(x="", y="", color = "Copy Number")+
    scale_x_continuous(expand = c(0.01, 0.01))
  
  pdf(paste0(o.dir, "qDNAseq_plots_updated/", sample.tmp, "_WGS_Scatterplot.pdf"), width = 5.5, height = 4)
  print(plotbulk)
  dev.off()
  
  
}
