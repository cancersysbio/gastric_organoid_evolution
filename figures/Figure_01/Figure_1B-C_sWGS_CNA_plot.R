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

sessionInfo()

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
wgs.files <- list.files("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure1/CNA_smooth_50kb_shallow", pattern = "CNA_smooth.txt", recursive = F, full.names = T, all.files = T)
all.wgs.files <- wgs.files


#####################################################################################
# FIGURE 1B - D1C1 TIMECOURSE
#####################################################################################
# set output directory
o.dir <- "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure1/PLOTS/Figure_1B_CNV/"
dir.create(paste0(o.dir, "qDNAseq_plots_updated"))
all.wgs.files
sample.ids <- c("D1")  
# add a sample id column
# read the data in
tables <- all.wgs.files[grep(sample.ids, all.wgs.files)]
tables <- tables[-grep("C2|C3|EARLY|LATE|Parent|P22|P26|87d|115d|582d", tables)]

tables
timepoints <- str_split_fixed(basename(tables), "_", 4)[,2]
timepoints
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
complete.wgs.data[complete.wgs.data$timepoint == "P3", "timepoint"] <- "WT"
complete.wgs.data[complete.wgs.data$timepoint == "55d", "timepoint"] <- "55 days"
complete.wgs.data[complete.wgs.data$timepoint == "101d", "timepoint"] <- "101 days"
complete.wgs.data[complete.wgs.data$timepoint == "190d", "timepoint"] <- "190 days"
complete.wgs.data[complete.wgs.data$timepoint == "260d", "timepoint"] <- "260 days"
complete.wgs.data[complete.wgs.data$timepoint == "344d", "timepoint"] <- "344 days"
complete.wgs.data[complete.wgs.data$timepoint == "442d", "timepoint"] <- "442 days"
complete.wgs.data[complete.wgs.data$timepoint == "608d", "timepoint"] <- "608 days"
complete.wgs.data$timepoint <- factor(complete.wgs.data$timepoint, levels = c("WT", "55 days", "101 days", "190 days", "260 days", "344 days", "442 days", "608 days"))
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
        axis.title = element_text(colour = "black", size = 24, face = "bold" ),
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

pdf(paste0(o.dir, "qDNAseq_plots_updated/", sample.ids, "_WGS_Scatterplot.pdf"), width = 7, height = 12)
print(plotbulk)
dev.off()

#####################################################################################
# FIGURE 1C - D1, D2 and D3 clones
#####################################################################################

o.dir <- "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure1/PLOTS/Figure_1C_CNV/"
dir.create(paste0(o.dir, "qDNAseq_plots_updated"))

sample.ids <- c("D1C1_608d", "D1C2_808d",  "D1C3_835d",
                "D2C1_588d", "D2C2_609d",  "D2C3_743d",
                "D3C1_718d", "D3C2_750d",  "D3C3_756d")  
# add a sample id column
# read the data in
tables <- all.wgs.files[grep(paste(sample.ids[1:9], collapse = "|"), all.wgs.files)]
clones <- str_split_fixed(basename(tables), "_", 6)[,1]
tables <- lapply(tables, read.table, stringsAsFactors = F, header = T)

complete.wgs.list <- list()
j <- 1
for (j in 1:length(tables)){
  
  table <- tables[[j]]
  colnames(table)[ncol(table)] <- "cnVal"
  table$clones <- clones[j]
  print(clones[j])
  
  # transform the data to get approx cnv states
  wgs.cnVal<- 2*2**table$cnVal
  wgs.cn.states <- wgs.cnVal
  
  # set up a dataframe with all information needed
  wgs.data <- data.frame(cn_state=wgs.cn.states, chr=table$chromosome, coordinates=c(1:50497), clones = table$clones)
  
  # Upperlimit the data to 8 copies
  wgs.data$cn_state[which(wgs.data$cn_state>8)] <- 8
  wgs.data$cn_state[which(wgs.data$cn_state < 0)] <- 0
  
  #This is needed to plot the chromosomes in the right order
  wgs.data$order <- c(1:50497)
  
  complete.wgs.list[[j]] <- wgs.data
}

complete.wgs.data <- bind_rows(complete.wgs.list)

donors <- c("D1", "D2", "D3")
donor <- donors[1]
for (donor in donors){
  
  donor.data <- complete.wgs.data[grep(donor, complete.wgs.data$clones),]
  donor.data <- donor.data[!is.na(donor.data$clones), ]
  donor.data$clones <- substr(donor.data$clones, 3,4)
  
  # create a ggplot for this
  plotbulk<- ggplot(donor.data, aes(coordinates, cn_state)) +
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
    facet_grid(clones ~ reorder(chr,order), scales="free", space="free", switch="both") + 
    ggExtra::removeGrid() + labs(x="Chromosomes", y="", color = "Copy Number")+
    scale_x_continuous(expand = c(0.01, 0.01))
  
  pdf(paste0(o.dir, "qDNAseq_plots_updated/", donor, "_WGS_Scatterplot.pdf"), width = 13.5, height = 4)
  print(plotbulk)
  dev.off()
  
}
