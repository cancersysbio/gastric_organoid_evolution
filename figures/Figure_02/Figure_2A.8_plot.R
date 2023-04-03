################################################################################################################################################
##                                                                                                                      
## CREATE MULTI COMPONENT FIGURE WITH SNVS, CNVS, SVS, SIGNATURES AND OTHER INFORMATION
##                                                                                                                      
##  Date: 10 SEPTEMBER 2021                                                                                                                    
##  
##  Author: Moritz Przybilla                                                                                                                    
##           
##                                                                                                                      
################################################################################################################################################

# clear workspace beforehand
rm(list = ls())

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("BiocManager", "reshape2", "ggrepel", "readr", "stringr", "tidyverse", "hdp", "sigfit", "BSgenome.Hsapiens.UCSC.hg38", "nrmisc", "lsa",
                      "lattice")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

#####################################################################################
# SET PATH
#####################################################################################


path <- "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_2/Figure_2A_my_reproduction"
setwd(file.path(path))



#####################################################################################
# READ IN MUTATION DATA
#####################################################################################

total_mutations <- read.table("Summary_table_Mutations.txt",header=TRUE)

head(total_mutations)


total_mutations <- subset(total_mutations, Sample != "D2C2_296d_MID")
total_mutations <- subset(total_mutations, Sample != "D3C1_296d_MID")

total_mutations$Patient <- factor(total_mutations$Patient, levels = c("D2", "D3"))
total_mutations$Annotation <- factor(total_mutations$Annotation, levels = c("EML", "Timecourse"))
total_mutations$Culture <- factor(total_mutations$Culture, levels = c("C1", "C2", "C3"))

total_mutations


#####################################################################################
# READ IN INDELS
#####################################################################################
indels.sample <- read.table("Summary_table_indels_old.txt",header=TRUE)
head(indels.sample)

indels.sample <- subset(indels.sample, Sample != "D2C2_296d_MID")
indels.sample <- subset(indels.sample, Sample != "D3C1_296d_MID")

indels.sample$Patient <- factor(indels.sample$Patient, levels = c("D2", "D3"))
indels.sample$Annotation <- factor(indels.sample$Annotation, levels = c("EML", "Timecourse"))
indels.sample$Culture <- factor(indels.sample$Culture, levels = c("C1", "C2", "C3"))


#####################################################################################
# READ IN WGII DATA
#####################################################################################

wgii.data.melt <- read.table("Summary_table_wGII.txt",header=TRUE)
head(wgii.data.melt)

wgii.data.melt$Patient <- factor(wgii.data.melt$Patient, levels = c("D2", "D3"))
wgii.data.melt$Annotation <- factor(wgii.data.melt$Annotation, levels = c("EML", "Timecourse"))
wgii.data.melt$Culture <- factor(wgii.data.melt$Culture, levels = c("C1", "C2", "C3"))


#####################################################################################
# READ IN LOH DATA
#####################################################################################

loh.data.melt <- read.table("Summary_table_LOH.txt",header=TRUE)
head(loh.data.melt)

loh.data.melt$Patient <- factor(loh.data.melt$Patient, levels = c("D2", "D3"))
loh.data.melt$Annotation <- factor(loh.data.melt$Annotation, levels = c("EML", "Timecourse"))
loh.data.melt$Culture <- factor(loh.data.melt$Culture, levels = c("C1", "C2", "C3"))


#####################################################################################
# READ IN SV CountS
#####################################################################################


sv.sample <- read.table("Summary_table_SV.txt",header=TRUE)
head(sv.sample)

sv.sample$Patient <- factor(sv.sample$Patient, levels = c("D2", "D3"))
sv.sample$Annotation <- factor(sv.sample$Annotation, levels = c("EML", "Timecourse"))
sv.sample$Culture <- factor(sv.sample$Culture, levels = c("C1", "C2", "C3"))


#####################################################################################
# CREATE THE PLOT
#####################################################################################

patients <- c("D2", "D3")
colfunc <- colorRampPalette(c("lightgrey", "darkblue"))

pdf("Figure_2A_no_Legend_15_10.pdf", width = 15, height = 10)
#pdf("test_mutations2.pdf", width = 15, height = 6)

#pdf("test_gene_loc3.pdf", width = 15, height = 6)


#pdf("test_full.pdf", width = 15, height = 40)


m <- matrix(c(1, 1, 2, 3, 3, 4, 5, 
              6, 6, 7, 8, 8, 9, 10,
              11, 11, 12, 13, 13, 14, 15,
              16, 16, 17, 18, 18, 19, 20, 
              21, 21, 22, 23, 23, 24, 25), nrow = 5, byrow = T)

par(oma = c(1, 3, 1, 1), mar=c(0.25, 1, 1, 1), mai = c(0.05, 0.05, 0.05, 0.25)) # bottom, left, top, right
layout(m)

#####################################################################################
#                       VISUALIZE SINGLE NUCLEOTIDE VARIANTS
#####################################################################################


i <- 1
for (i in 1:length(patients)){
  
  # which patient
  patient <- patients[i]
  print("patient")
  print(patient)
  
  # get dataframe subset
  df.subset <- total_mutations[total_mutations$Patient == patient, ]
  Annotation <- c("Timecourse", "EML")
  
  j <- 1
  
  for (j in 1:length(Annotation)){
    
    # which Annotation
    Annotation.of.interest <- Annotation[j]
    print(Annotation.of.interest)
    
    #
    df.Annotation.subset <- df.subset[df.subset$Annotation == Annotation.of.interest, ]
    Cultures <- unique(df.Annotation.subset$Culture)
    k <- 1
    
    for (k in 1:length(Cultures)){
      
      # which Culture
      Culture.of.interest <- Cultures[k]
      print("Culture.of.interest")
      print(Culture.of.interest)
      
      #
      final.subset <- df.Annotation.subset[df.Annotation.subset$Culture == Culture.of.interest, ]
      print ("final.subset")
      print (final.subset)
      
      if (Annotation.of.interest == "Timecourse"){
        
        final.subset$Timepoint <- factor(final.subset$Timepoint, levels = c("115d", "173d", "190d", "260d", "264d", "404d", "428d", "718d", "729d"))
        print(final.subset)
      } else if(Annotation.of.interest == "EML"){
        
        final.subset$Timepoint <- factor(final.subset$Timepoint, levels = c("MID", "LATE"))
        
      }
      
      if (Annotation.of.interest == "EML"){
        
        barplot(final.subset$Count, ylim = c(0, 14000), #main = Culture.of.interest,
                ylab = "SNV", axes = F,  axisnames=F, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, width=c(5, 5),
                col = c("darkblue", "darkmagenta"))
                #legend.text = final.subset$Timepoint) # Stacked bars (default)
        box()
        
      } else {
        
        barplot(final.subset$Count,  ylim = c(0, 14000), #main = Culture.of.interest, 
                ylab = "SNV", axes = F, axisnames=F, cex.lab=1.5, cex.axis=1.5, cex.main=1.5,
                col =  colfunc(5))
                #legend.text = final.subset$Timepoint)
        box()
        
      }
      
      # add axis only to the first plot
      if (i == 1 & j == 1 & k == 1){
        axis(2, font=1.5)
      }
    }
  }
}


#####################################################################################
#                               VISUALIZE INDELS
#####################################################################################


i <- 1
for (i in 1:length(patients)){
  
  # which patient
  patient <- patients[i]
  print(patient)
  
  # get dataframe subset
  df.subset <- indels.sample[indels.sample$Patient == patient, ]
  Annotation <- c("Timecourse", "EML")
  j <- 1
  
  for (j in 1:length(Annotation)){
    
    # which Annotation
    Annotation.of.interest <- Annotation[j]
    print(Annotation.of.interest)
    
    #
    df.Annotation.subset <- df.subset[df.subset$Annotation == Annotation.of.interest, ]
    Cultures <- unique(df.Annotation.subset$Culture)
    k <- 1
    
    for (k in 1:length(Cultures)){
      
      # which Culture
      Culture.of.interest <- Cultures[k]
      print(Culture.of.interest)
      
      #
      final.subset <- df.Annotation.subset[df.Annotation.subset$Culture == Culture.of.interest, ]
      
      if (Annotation.of.interest == "Timecourse"){
        
        final.subset$Timepoint <- factor(final.subset$Timepoint, levels = c("115d", "173d", "190d", "260d", "264d", "404d", "428d", "718d", "729d"))
        
      }
      
      # reshape it
      matrix <-acast(final.subset, final.subset$Sample ~ final.subset$Type, value.var = "Count")
      
      # # order according to mutations
      # final.subset <- final.subset[order(final.subset$total_muts),]
      
      if (Annotation.of.interest == "EML"){
        
        matrix <- matrix[c(grep("MID", rownames(matrix)),grep("LATE", rownames(matrix))),]
        
        barplot(t(matrix),
                ylab = "Indel",axes = F,  axisnames=F, ylim = c(0, 800), 
                #ylab = "Indel",axes = F,  axisnames=F, ylim = c(0, 1600), 
                
                col = c("darkorange", "purple"), cex.lab=1.5, cex.axis=1.5,
                #legend.text = colnames(matrix),
                beside = FALSE) # Stacked bars (default)
        
        
        box()
        
      } else {
        
        barplot(t(matrix), ylim = c(0, 800), 
        #barplot(t(matrix), ylim = c(0, 1600), 
                ylab = "CNV", axes = F,  axisnames=F,
                col = c("darkorange", "purple"), cex.lab=1.5, cex.axis=1.5,
                #legend.text = colnames(matrix),
                beside = FALSE) # Stacked bars (default)
        box()
        
      }
      
      # add axis only to the first plot
      if (i == 1 & j == 1 & k == 1){
        axis(2, font=1.5)
      }
    }
  }
}

#####################################################################################
#                            VISUALIZE WGII
#####################################################################################

i <- 1
for (i in 1:length(patients)){
  
  # which patient
  patient <- patients[i]
  print(patient)
  
  # get dataframe subset
  df.subset <- wgii.data.melt[wgii.data.melt$Patient == patient, ]
  Annotation <- c("Timecourse", "EML")
  j <- 1
  
  for (j in 1:length(Annotation)){
    
    # which Annotation
    Annotation.of.interest <- Annotation[j]
    print(Annotation.of.interest)
    
    #
    df.Annotation.subset <- df.subset[df.subset$Annotation == Annotation.of.interest, ]
    Cultures <- unique(df.Annotation.subset$Culture)
    k <- 1
    
    for (k in 1:length(Cultures)){
      
      # which Culture
      Culture.of.interest <- Cultures[k]
      print(Culture.of.interest)
      
      final.subset <- df.Annotation.subset[df.Annotation.subset$Culture == Culture.of.interest, ]
      
      if (Annotation.of.interest == "Timecourse"){
        
        final.subset$Timepoint <- factor(final.subset$Timepoint, levels = c("115d", "173d", "190d", "260d", "264d", "404d", "428d", "718d", "729d"))
        
      } else if(Annotation.of.interest == "EML"){
        
        final.subset$Timepoint <- factor(final.subset$Timepoint, levels = c("MID", "LATE"))
        
      }
      
      if (Annotation.of.interest == "EML"){
        
        barplot(final.subset$Fraction, ylim = c(0, 0.2), #main = Culture.of.interest,
                ylab = "wGII", axes = F,  axisnames=F, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, width=c(5, 5),
                col = c("darkblue", "darkmagenta"))
                #legend.text = final.subset$Timepoint) # Stacked bars (default)
        box()
        
      } else {
        
        barplot(final.subset$Fraction,  ylim = c(0, 0.2), #main = Culture.of.interest, 
                ylab = "wGII", axes = F, axisnames=F, cex.lab=1.5, cex.axis=1.5, cex.main=1.5,
                col =  colfunc(5))
                #legend.text = final.subset$Timepoint)
        box()
        
      }
      
      # add axis only to the first plot
      if (i == 1 & j == 1 & k == 1){
        axis(2, font=1.5)
      }
    }
  }
}

?barplot


#####################################################################################
#                            VISUALIZE LOH
#####################################################################################

i <- 1
for (i in 1:length(patients)){
  
  # which patient
  patient <- patients[i]
  print(patient)
  
  # get dataframe subset
  df.subset <- loh.data.melt[loh.data.melt$Patient == patient, ]
  Annotation <- c("Timecourse", "EML")
  j <- 1
  
  for (j in 1:length(Annotation)){
    
    # which Annotation
    Annotation.of.interest <- Annotation[j]
    print(Annotation.of.interest)
    
    #
    df.Annotation.subset <- df.subset[df.subset$Annotation == Annotation.of.interest, ]
    Cultures <- unique(df.Annotation.subset$Culture)
    k <- 1
    
    for (k in 1:length(Cultures)){
      
      # which Culture
      Culture.of.interest <- Cultures[k]
      print(Culture.of.interest)
      
      final.subset <- df.Annotation.subset[df.Annotation.subset$Culture == Culture.of.interest, ]
      
      if (Annotation.of.interest == "Timecourse"){
        
        final.subset$Timepoint <- factor(final.subset$Timepoint, levels = c("115d", "173d", "190d", "260d", "264d", "404d", "428d", "718d", "729d"))
        
      } else if(Annotation.of.interest == "EML"){
        
        final.subset$Timepoint <- factor(final.subset$Timepoint, levels = c("MID", "LATE"))
        
      }
      
      if (Annotation.of.interest == "EML"){
        
        barplot(final.subset$Fraction, ylim = c(0, 0.2), #main = Culture.of.interest,
                ylab = "LOH", axes = F,  axisnames=F, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, width=c(5, 5),
                col = c("darkblue", "darkmagenta"))
                #legend.text = final.subset$Timepoint) # Stacked bars (default)
        box()
        
      } else {
        
        barplot(final.subset$Fraction,  ylim = c(0, 0.2), #main = Culture.of.interest, 
                ylab = "LOH", axes = F, axisnames=F, cex.lab=1.5, cex.axis=1.5, cex.main=1.5,
                col =  colfunc(5))
                #legend.text = final.subset$Timepoint)
        box()
        
      }
      
      # add axis only to the first plot
      if (i == 1 & j == 1 & k == 1){
        axis(2, font=1.5)
      }
    }
  }
}

#####################################################################################
#                       VISUALIZE STRUCTURAL VARIANTS
#####################################################################################

i <- 1
for (i in 1:length(patients)){
  
  # which patient
  patient <- patients[i]
  print(patient)
  
  # get dataframe subset
  df.subset <- sv.sample[sv.sample$Patient == patient, ]
  Annotation <- c("Timecourse", "EML")
  j <- 1
  
  for (j in 1:length(Annotation)){
    
    # which Annotation
    Annotation.of.interest <- Annotation[j]
    print(Annotation.of.interest)
    
    #
    df.Annotation.subset <- df.subset[df.subset$Annotation == Annotation.of.interest, ]
    Cultures <- unique(df.Annotation.subset$Culture)
    k <- 1
    
    for (k in 1:length(Cultures)){
      
      # which Culture
      Culture.of.interest <- Cultures[k]
      print(Culture.of.interest)
      
      #
      final.subset <- df.Annotation.subset[df.Annotation.subset$Culture == Culture.of.interest, ]
      
      if (Annotation.of.interest == "Timecourse"){
        
        final.subset$Timepoint <- factor(final.subset$Timepoint, levels = c("115d", "173d", "190d", "260d", "264d", "404d", "428d", "718d", "729d"))
        
      }
      
      # reshape it
      matrix <-acast(final.subset, final.subset$Sample ~ final.subset$Type, value.var = "Count")
      matrix[is.na(matrix)] <- "0"
      
      # # order according to mutations
      # final.subset <- final.subset[order(final.subset$total_muts),]
      
      if (Annotation.of.interest == "EML"){
        
        matrix <- matrix[c(grep("MID", rownames(matrix)),grep("LATE", rownames(matrix))),]
        
        barplot(t(matrix),
                ylab = "SV",axes = F,  axisnames=F, ylim = c(0, 100), 
                col = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"), cex.lab=1.5, cex.axis=1.5,
                #legend.text = colnames(matrix),
                beside = FALSE) # Stacked bars (default)
        
        
        box()
        
      } else {
        
        barplot(t(matrix), ylim = c(0, 100), 
                ylab = "SV", axes = F,  axisnames=F,
                col = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"), cex.lab=1.5, cex.axis=1.5,
                #legend.text = colnames(matrix),
                beside = FALSE) # Stacked bars (default)
        
        
        box()
        
      }
      
      # add axis only to the first plot
      if (i == 1 & j == 1 & k == 1){
        axis(2, font=1.5)
      }
    }
  }
}
dev.off()

