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
wgs.files <- list.files("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_2/INPUT/CNA_smooth_50kb_deep", pattern = "CNA_smooth.txt", recursive = F, full.names = T, all.files = T)
all.wgs.files <- wgs.files

o.dir <- "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_2/PLOTS/Figure_2D_CNV/"
dir.create(paste0(o.dir, "qDNAseq_plots_updated"))


#####################################################################################
# FIGURE 2B D3C1 TIMECOURSE DEEP SEQUENCING
#####################################################################################


sample.ids <- c("D3C1")  
# add a sample id column
# read the data in
tables <- all.wgs.files[grep(sample.ids, all.wgs.files)]
tables <- tables[-grep("D3C1_296d_MID", tables)]
timepoints <- str_split_fixed(basename(tables), "_", 6)[,2]
tables <- lapply(tables, read.table, stringsAsFactors = F, header = T)
complete.wgs.list <- list()

j <- 1
for (j in 1:length(tables)){
  
  table <- tables[[j]]
  colnames(table)[ncol(table)] <- "cnVal"
  table$timepoint <- timepoints[j]
  print(timepoints[j])
  
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
  
  complete.wgs.list[[j]] <- wgs.data
}

complete.wgs.data <- bind_rows(complete.wgs.list)

# replace weeks with days
complete.wgs.data[complete.wgs.data$timepoint == "115d", "timepoint"] <- "115 days"
complete.wgs.data[complete.wgs.data$timepoint == "173d", "timepoint"] <- "173 days"
complete.wgs.data[complete.wgs.data$timepoint == "264d", "timepoint"] <- "264 days"
complete.wgs.data[complete.wgs.data$timepoint == "404d", "timepoint"] <- "404 days"
complete.wgs.data[complete.wgs.data$timepoint == "718d", "timepoint"] <- "718 days"
complete.wgs.data$timepoint <- factor(complete.wgs.data$timepoint, levels = c("115 days", "173 days", "264 days", "404 days", "718 days"))
complete.wgs.data <- complete.wgs.data[!is.na(complete.wgs.data$timepoint), ]

# create a ggplot for this
plotbulk<- ggplot(complete.wgs.data, aes(coordinates, cn_state)) +
  geom_point(aes(colour = cn_state), size = 0.01) +
  scale_colour_gradient2(low = "midnightblue", high = "deeppink4", mid ="#F8F8F8", midpoint = 2, limits=c(0,4)) +
  theme_bw() +
  scale_y_continuous(breaks = c(0, 2, 4), limits = c(0,4), position = "right") +
  theme(panel.grid.major = element_blank(),
        strip.text = element_text(face="bold", size=20, colour = "black",),
        strip.background = element_rect(fill="white", colour="black", size=1), 
        strip.text.x = element_blank(),
        axis.text = element_text(colour = "black", size = 14, face = "bold" ),
        axis.title = element_text(colour = "black", size = 16, face = "bold" ),
        plot.title = element_text(colour = "black", size = 12, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 9, face = "bold",),
        legend.text = element_text(colour="black", size=7, face="bold"),
        legend.position = "None",
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.spacing.x=unit(0, "lines"), panel.border = element_rect(linetype =3))+
  facet_grid(timepoint ~ reorder(chr,order), scales="free", space="free", switch="both") + 
  ggExtra::removeGrid() + labs(x="Chromosomes", y="", color = "Copy Number")+
  scale_x_continuous(expand = c(0.01, 0.01))

pdf(paste0(o.dir, "qDNAseq_plots_updated/", sample.ids, "_WGS_Scatterplot_10x10.pdf"), width = 10, height = 10)
print(plotbulk )
dev.off()



#####################################################################################
# SUPPLEMENTARY FIGURE 2B D2C2 TIMECOURSE DEEP SEQUENCING
#####################################################################################

o.dir <- "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_2/PLOTS/Figure_ED_4C_CNV/"
dir.create(paste0(o.dir, "qDNAseq_plots_updated"))

sample.ids <- c("D2C2")  
# add a sample id column
# read the data in
tables <- all.wgs.files[grep(sample.ids, all.wgs.files)]
tables <- tables[-grep("D2C2_296d_MID|R2T", tables)]
timepoints <- str_split_fixed(basename(tables), "_", 6)[,2]
tables <- lapply(tables, read.table, stringsAsFactors = F, header = T)
complete.wgs.list <- list()

j <- 1
for (j in 1:length(tables)){
  
  table <- tables[[j]]
  colnames(table)[ncol(table)] <- "cnVal"
  table$timepoint <- timepoints[j]
  print(timepoints[j])
  
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
  
  complete.wgs.list[[j]] <- wgs.data
}



complete.wgs.data <- bind_rows(complete.wgs.list)

# replace weeks with days
complete.wgs.data[complete.wgs.data$timepoint == "115d", "timepoint"] <- "115 days"
complete.wgs.data[complete.wgs.data$timepoint == "190d", "timepoint"] <- "190 days"
complete.wgs.data[complete.wgs.data$timepoint == "260d", "timepoint"] <- "260 days"
complete.wgs.data[complete.wgs.data$timepoint == "428d", "timepoint"] <- "428 days"
complete.wgs.data[complete.wgs.data$timepoint == "729d", "timepoint"] <- "729 days"
complete.wgs.data$timepoint <- factor(complete.wgs.data$timepoint, levels = c("115 days", "190 days", "260 days", "428 days", "729 days"))
complete.wgs.data <- complete.wgs.data[!is.na(complete.wgs.data$timepoint), ]

# create a ggplot for this
plotbulk<- ggplot(complete.wgs.data, aes(coordinates, cn_state)) +
  geom_point(aes(colour = cn_state), size = 0.01) +
  scale_colour_gradient2(low = "midnightblue", high = "deeppink4", mid ="#F8F8F8", midpoint = 2, limits=c(0,4)) +
  theme_bw() +
  scale_y_continuous(breaks = c(0, 2, 4), limits = c(0,4), position = "right") +
  theme(panel.grid.major = element_blank(),
        strip.text = element_text(face="bold", size=20, colour = "black",),
        strip.background = element_rect(fill="white", colour="black", size=1), 
        strip.text.x = element_blank(),
        axis.text = element_text(colour = "black", size = 14, face = "bold" ),
        axis.title = element_text(colour = "black", size = 16, face = "bold" ),
        plot.title = element_text(colour = "black", size = 12, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 9, face = "bold",),
        legend.text = element_text(colour="black", size=7, face="bold"),
        legend.position = "None",
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.spacing.x=unit(0, "lines"), panel.border = element_rect(linetype =3))+
  facet_grid(timepoint ~ reorder(chr,order), scales="free", space="free", switch="both") + 
  ggExtra::removeGrid() + labs(x="Chromosomes", y="", color = "Copy Number")+
  scale_x_continuous(expand = c(0.01, 0.01))

pdf(paste0(o.dir, "qDNAseq_plots_updated/", sample.ids, "_WGS_Scatterplot_10x10.pdf"), width = 10, height = 10)
print(plotbulk)
dev.off()

#####################################################################################
# FIGURE 2 FHIT LOCUS CNV PLOTS 
#####################################################################################

wgs.files <- list.files("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_2/INPUT/CNA_smooth_1kb_deep", recursive = T, full.names = T, all.files = T)
all.wgs.files <- wgs.files

o.dir <- "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_2/PLOTS/Figure_2E_CNV/"
dir.create(paste0(o.dir, "qDNAseq_plots_updated"))

# add a sample id column
# read the data in
sample.ids <- c("D3C1")  
tables <- all.wgs.files[grep(sample.ids, all.wgs.files)]
tables <- tables[-grep("D3C1_296d_MID|chrarm", tables)]
sample.ids <- str_split_fixed(basename(tables), "_", 6)[,2]
timepoints <- str_split_fixed(basename(tables), "_", 6)[,2]
complete.wgs.list <- list()



# add a sample id column
# read the data in
tables <- lapply(tables, read.table, stringsAsFactors = F, header = T)

head(tables)
head(table)

j <- 1
for (j in 1:length(tables)){
  table <- tables[[j]]
  colnames(table)[ncol(table)] <- "cnVal"
  table$timepoint <- timepoints[j]
  print(timepoints[j])
  
  # reduce to chromosome 3 60MB - 60.1MB
  table <- table[table$chromosome == 3 & table$start >= 60000000 & table$end <= 61000000,]
  
  # transform the data to get approx cnv states
  wgs.cnVal<- 2*2**table$cnVal
  wgs.cn.states <- wgs.cnVal
  
  # set up a dataframe with all information needed
  wgs.data <- data.frame(cn_state=wgs.cn.states, chr=table$chromosome, coordinates=c(1:length(wgs.cn.states)), timepoint = table$timepoint)
  
  # Upperlimit the data to 8 copies
  wgs.data$cn_state[which(wgs.data$cn_state>8)] <- 8
  wgs.data$cn_state[which(wgs.data$cn_state < 0)] <- 0
  
  #This is needed to plot the chromosomes in the right order
  wgs.data$order <- c(1:length(wgs.cn.states))
  
  complete.wgs.list[[j]] <- wgs.data
}

complete.wgs.list

complete.wgs.data <- bind_rows(complete.wgs.list)
hist(complete.wgs.data$cn_state, breaks = 1000)

i <- 5
for (i in 1:length(sample.ids)){
  
  sample.tmp <- sample.ids[i]
  print(sample.tmp)
  sample.wgs.data <- complete.wgs.data[complete.wgs.data$timepoint == sample.tmp,]
  print(nrow(sample.wgs.data))
  
  # create a ggplot for this
  plotbulk<- ggplot(sample.wgs.data, aes(coordinates, cn_state)) +
    geom_point(aes(colour = cn_state), size = 0.3) +
    scale_colour_gradient2(low = "midnightblue", high = "deeppink4", mid ="#F8F8F8", midpoint = 2, limits=c(0,4),oob = scales::squish) +
    theme_bw() +
    scale_y_continuous(breaks = c(0, 2, 4), limits = c(0,4), position = "left") +   ### MAYBE CHANGE
    theme(panel.grid.major = element_blank(),
          strip.text = element_text(face="bold", size=20, colour = "black",),
          #strip.background = element_rect(fill="white", colour="black", size=1), 
          strip.text.x = element_blank(),
          axis.line = element_line(size = 0.5, colour = "black"),
          axis.ticks.y = element_line(size = 0.5, colour = "black"),
          axis.text = element_text(colour = "black", size = 14, face = "bold" ),
          axis.title = element_text(colour = "black", size = 16, face = "bold" ),
          plot.title = element_text(colour = "black", size = 12, face = "bold", hjust = 0.5),
          legend.title = element_text(color = "black", size = 9, face = "bold",),
          legend.text = element_text(colour="black", size=7, face="bold"),
          legend.position = "None",
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "white", fill=NA, size=1),
          #panel.spacing.x=unit(0, "lines"),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())+
    facet_grid( ~ reorder(chr,order), scales="free", space="free", switch="both") + 
    ggExtra::removeGrid() + labs(x="", y="", color = "Copy Number")+
    scale_x_continuous(expand = c(0.01, 0.01))
  
  pdf(paste0(o.dir, "qDNAseq_plots_updated/", sample.ids[i], "_WGS_FHIT_Scatterplot.pdf"), width = 5.5, height = 1.1)
  print(plotbulk)
  dev.off()
  
}



#####################################################################################
# SUPPLEMENTARY FIGURE 2 PYRGO LOCUS CNV PLOTS 
#####################################################################################

wgs.files <- list.files("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_2/INPUT/CNA_smooth_1kb_deep", recursive = T, full.names = T, all.files = T)
all.wgs.files <- wgs.files

o.dir <- "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_2/PLOTS/Figure_ED_4B_CNV/"

dir.create(paste0(o.dir, "qDNAseq_plots_updated"))

# add a sample id column
# read the data in
sample.ids <- c("D2C2")  
tables <- all.wgs.files[grep(sample.ids, all.wgs.files)]
tables <- tables[-grep("296d_MID|R2T|chrarm", tables)]
sample.ids <- str_split_fixed(basename(tables), "_", 6)[,2]
timepoints <- str_split_fixed(basename(tables), "_", 6)[,2]
complete.wgs.list <- list()

tables
timepoints
# add a sample id column
# read the data in
tables <- lapply(tables, read.table, stringsAsFactors = F, header = T)

j <- 1
for (j in 1:length(tables)){
  
  table <- tables[[j]]
  colnames(table)[ncol(table)] <- "cnVal"
  table$timepoint <- timepoints[j]
  print(timepoints[j])
  
  # reduce to chromosome 3 60MB - 60.1MB
  table <- table[table$chromosome == 3 & table$start >= 86000000 & table$end <= 88000000,]
  
  # transform the data to get approx cnv states
  wgs.cnVal<- 2*2**table$cnVal
  wgs.cn.states <- wgs.cnVal
  
  # set up a dataframe with all information needed
  wgs.data <- data.frame(cn_state=wgs.cn.states, chr=table$chromosome, coordinates=c(1:length(wgs.cn.states)), timepoint = table$timepoint)
  
  # Upperlimit the data to 8 copies
  wgs.data$cn_state[which(wgs.data$cn_state>8)] <- 8
  wgs.data$cn_state[which(wgs.data$cn_state < 0)] <- 0
  
  #This is needed to plot the chromosomes in the right order
  wgs.data$order <- c(1:length(wgs.cn.states))
  
  complete.wgs.list[[j]] <- wgs.data
}

complete.wgs.data <- bind_rows(complete.wgs.list)
hist(complete.wgs.data$cn_state, breaks = 1000)

i <- 5
for (i in 1:length(sample.ids)){
  
  sample.tmp <- sample.ids[i]
  print(sample.tmp)
  sample.wgs.data <- complete.wgs.data[complete.wgs.data$timepoint == sample.tmp,]
  print(nrow(sample.wgs.data))
  
  # create a ggplot for this
  plotbulk<- ggplot(sample.wgs.data, aes(coordinates, cn_state)) +
    geom_point(aes(colour = cn_state), size = 0.3) +
    scale_colour_gradient2(low = "midnightblue", high = "deeppink4", mid ="#F8F8F8", midpoint = 2, limits=c(0,4),oob = scales::squish) +
    theme_bw() +
    scale_y_continuous(breaks = c(0, 2, 4, 6, 8), limits = c(0,8), position = "left") +   ### MAYBE CHANGE
    theme(panel.grid.major = element_blank(),
          strip.text = element_text(face="bold", size=20, colour = "black",),
          #strip.background = element_rect(fill="white", colour="black", size=1), 
          strip.text.x = element_blank(),
          axis.line = element_line(size = 0.5, colour = "black"),
          axis.ticks.y = element_line(size = 0.5, colour = "black"),
          axis.text = element_text(colour = "black", size = 14, face = "bold" ),
          axis.title = element_text(colour = "black", size = 16, face = "bold" ),
          plot.title = element_text(colour = "black", size = 12, face = "bold", hjust = 0.5),
          legend.title = element_text(color = "black", size = 9, face = "bold",),
          legend.text = element_text(colour="black", size=7, face="bold"),
          legend.position = "None",
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "white", fill=NA, size=1),
          #panel.spacing.x=unit(0, "lines"),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())+
    facet_grid( ~ reorder(chr,order), scales="free", space="free", switch="both") + 
    ggExtra::removeGrid() + labs(x="", y="", color = "Copy Number")+
    scale_x_continuous(expand = c(0.01, 0.01))
  
  pdf(paste0(o.dir, "qDNAseq_plots_updated/", sample.ids[i], "_WGS_PYRGO_Scatterplot_8.pdf"), width = 5.5, height = 1.1)
  print(plotbulk)
  dev.off()
  
}

