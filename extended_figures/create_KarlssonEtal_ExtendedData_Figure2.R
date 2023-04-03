################################################################################################################################################
##                                                                                                                      
##  CREATE MULTI COMPONENT SUPPLEMENTARY FIGURE WITH SNVS, CNVS, SVS, SIGNATURES AND OTHER INFORMATION
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
list.of.packages <- c("BiocManager", "reshape2", "ggrepel", "readr", "stringr", "tidyverse", "lattice", "data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

#####################################################################################
# READ INPUT DATA
#####################################################################################

cnv.melt.summary <- fread("/Users/mp34/stanford/Data/Correct_files/Summary_table_CNA.txt", header = T, sep = "\t", data.table = F)
total_mutations <- fread("/Users/mp34/stanford/Data/Correct_files/Summary_table_mutations.txt", header = T, sep = "\t", data.table = F)
wgii.data.melt <- fread("/Users/mp34/stanford/Data/Correct_files/Summary_table_wGII.txt", header = T, sep = "\t", data.table = F)
loh.data.melt <- fread("/Users/mp34/stanford/Data/Correct_files/Summary_table_LOH.txt", header = T, sep = "\t", data.table = F)
sample_signatures.melt <- fread("/Users/mp34/stanford/Data/Correct_files/Summary_table_signatures.txt", header = T, sep = "\t", data.table = F)
cosmic.sample <- fread("/Users/mp34/stanford/Data/Correct_files/Summary_table_cosmic.txt", header = T, sep = "\t", data.table = F)
gene.mutation.sample <- fread("/Users/mp34/stanford/Data/Correct_files/Summary_table_location.txt", header = T, sep = "\t", data.table = F)
kataegis.sample <- fread("/Users/mp34/stanford/Data/Correct_files/Summary_table_kataegis.txt", header = T, sep = "\t", data.table = F)
indels.sample <- fread("/Users/mp34/stanford/Data/Correct_files/Summary_table_indels_old.txt", header = T, sep = "\t", data.table = F)
sv.sample <- fread("/Users/mp34/stanford/Data/Correct_files/Summary_table_SV.txt", header = T, sep = "\t", data.table = F)
sv.prop <- fread("/Users/mp34/stanford/Data/Correct_files/Summary_table_SV_prop.txt", header = T, sep = "\t", data.table = F)

#####################################################################################
# CREATE THE PLOT
#####################################################################################

patients <- c("D2", "D3")
colfunc <- colorRampPalette(c("lightgrey", "darkblue"))

pdf("/Users/mp34/stanford/KarlssonEtal_ExtendedFigure.pdf", width = 15, height = 10)

m <- matrix(c(1, 1, 2, 3, 3, 4, 5, 
              6, 6, 7, 8, 8, 9, 10,
              11, 11, 12, 13, 13, 14, 15,
              16, 16, 17, 18, 18, 19, 20, 
              21, 21, 22, 23, 23, 24, 25,
              26, 26, 27, 28, 28, 29, 30,
              31, 31, 32, 33, 33, 34, 35,
              36, 36, 37, 38, 38, 39, 40,
              41, 41, 42, 43, 43, 44, 45,
              46, 46, 47, 48, 48, 49, 50,
              51, 51, 52, 53, 53, 54, 55), nrow = 11, byrow = T)

par(oma = c(1, 3, 1, 1), mar=c(0.25, 1, 1, 1), mai = c(0.05, 0.05, 0.05, 0.25)) # bottom, left, top, right
layout(m)

#####################################################################################
#                       VISUALIZE SINGLE NUCLEOTIDE VARIANTS
#####################################################################################

i <- 1
for (i in 1:length(patients)){
  
  # which patient
  patient <- patients[i]
  print(patient)
  
  # get dataframe subset
  df.subset <- total_mutations[total_mutations$Patient == patient, ]
  anno <- c("Timecourse", "EML")
  
  j <- 1
  
  for (j in 1:length(anno)){
    
    # which anno
    anno.of.interest <- anno[j]
    print(anno.of.interest)
    
    #
    df.anno.subset <- df.subset[df.subset$Annotation == anno.of.interest, ]
    clones <- unique(df.anno.subset$Culture)
    k <- 1
    
    for (k in 1:length(clones)){
      
      # which clone
      clone.of.interest <- clones[k]
      print(clone.of.interest)
      
      #
      final.subset <- df.anno.subset[df.anno.subset$Culture == clone.of.interest, ]
      
      if (anno.of.interest == "Timecourse"){
        
        final.subset$Timepoint <- factor(final.subset$Timepoint, levels = c("101d", "115d", "173d", "190d", "260d", "264d", "404d", "428d", "609d", "705d", "718d", "729d"))
        
      } else if(anno.of.interest == "EML"){
        
        final.subset$Timepoint <- factor(final.subset$Timepoint, levels = c("MID", "LATE"))
        
      }
      
      if (anno.of.interest == "EML"){
        
        barplot(final.subset$Count, ylim = c(0, 12500), main = clone.of.interest,
                ylab = "SNV", axes = F,  axisnames=F, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, width=c(5, 5),
                col = c("darkblue", "darkmagenta"),
                legend.text = final.subset$Timepoint) # Stacked bars (default)
        box()
        
      } else {
        
        barplot(final.subset$Count,  ylim = c(0, 12500), main = clone.of.interest, axes = F,
                ylab = "SNV", axisnames=F, cex.lab=1.5, cex.axis=1.5, cex.main=1.5,
                col =  colfunc(5),
                legend.text = final.subset$Timepoint)
        box()
        
      }
      
      # add axis only to the first plot
      if (i == 1 & j == 1 & k == 1){
        axis(2, font=1.5)
      }
    }
  }
}

# ####################################################################################
#                             VISUALIZE SIGNATURES
# ####################################################################################

i <- 1
for (i in 1:length(patients)){
  
  # which patient
  patient <- patients[i]
  print(patient)
  
  # get dataframe subset
  df.subset <- sample_signatures.melt[sample_signatures.melt$Patient == patient, ]
  anno <- c("Timecourse", "EML")
  j <- 1
  
  for (j in 1:length(anno)){
    
    # which anno
    anno.of.interest <- anno[j]
    print(anno.of.interest)
    
    #
    df.anno.subset <- df.subset[df.subset$Annotation == anno.of.interest, ]
    clones <- unique(df.anno.subset$Culture)
    k <- 1
    
    for (k in 1:length(clones)){
      
      # which clone
      clone.of.interest <- clones[k]
      print(clone.of.interest)
      
      #
      final.subset <- df.anno.subset[df.anno.subset$Culture == clone.of.interest, ]
      
      if (anno.of.interest == "Timecourse"){
        
        final.subset$Timepoint <- factor(final.subset$Timepoint, levels = c("101d", "115d", "173d", "190d", "260d", "264d", "404d", "428d", "609d", "705d", "718d", "729d"))
        
      } else if(anno.of.interest == "EML"){
        
        final.subset$Timepoint <- factor(final.subset$Timepoint, levels = c("MID", "LATE"))
        
      }
      
      # reshape it
      matrix <- acast(final.subset, final.subset$Sample ~ final.subset$Signature, value.var = "Fraction")
      matrix <- matrix[, c("SBS1", "SBS5", "SBS17a", "SBS17b", "SBS40")]
      
      # # order according to mutations
      # final.subset <- final.subset[order(final.subset$total_muts),]
      
      if (anno.of.interest == "EML"){
        
        matrix <- matrix[c(grep("MID", rownames(matrix)),grep("LATE", rownames(matrix))),]
        
        barplot(t(matrix),
                ylab = "Signatures",axes = F,  axisnames=F, ylim = c(0, 1),
                col=c(RColorBrewer::brewer.pal(12, "Paired"),"magenta","firebrick"), cex.lab=1.5, cex.axis=1.5,
                legend.text = colnames(matrix),
                beside = FALSE) # Stacked bars (default)
        
        
        box()
        
      } else {
        
        barplot(t(matrix), ylim = c(0, 1),
                ylab = "Signatures", axes = F,  axisnames=F,
                col=c(RColorBrewer::brewer.pal(12, "Paired"),"magenta","firebrick"), cex.lab=1.5, cex.axis=1.5,
                legend.text = colnames(matrix),
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

# #####################################################################################
# #               VISUALIZE SINGLE NUCLEOTIDE VARIANTS IN COSMIC
# #####################################################################################

i <- 1
for (i in 1:length(patients)){
  
  # which patient
  patient <- patients[i]
  print(patient)
  
  # get dataframe subset
  df.subset <- cosmic.sample[cosmic.sample$Patient == patient, ]
  anno <- c("Timecourse", "EML")
  j <- 1
  
  for (j in 1:length(anno)){
    
    # which anno
    anno.of.interest <- anno[j]
    print(anno.of.interest)
    
    #
    df.anno.subset <- df.subset[df.subset$Annotation == anno.of.interest, ]
    clones <- unique(df.anno.subset$Culture)
    k <- 1
    
    for (k in 1:length(clones)){
      
      # which clone
      clone.of.interest <- clones[k]
      print(clone.of.interest)
      
      #
      final.subset <- df.anno.subset[df.anno.subset$Culture == clone.of.interest, ]
      
      if (anno.of.interest == "Timecourse"){
        
        final.subset$Timepoint <- factor(final.subset$Timepoint, levels = c("101d", "115d", "173d", "190d", "260d", "264d", "404d", "428d", "609d", "705d", "718d", "729d"))
        
      } else if(anno.of.interest == "EML"){
        
        final.subset$Timepoint <- factor(final.subset$Timepoint, levels = c("MID", "LATE"))
        
      }
      
      if (anno.of.interest == "EML"){
        
        barplot(final.subset$Count, ylim = c(0, 400), main = clone.of.interest,
                ylab = "SNV", axes = F,  axisnames=F, cex.lab=1.5, cex.axis=1.5, cex.main=1.5,
                col = c("darkblue", "darkmagenta"),
                legend.text = final.subset$Timepoint) # Stacked bars (default)
        box()
        
      } else {
        
        barplot(final.subset$Count,  ylim = c(0, 400), main = clone.of.interest, axes = F,
                ylab = "SNV", axisnames=F, cex.lab=1.5, cex.axis=1.5, cex.main=1.5,
                col =   colfunc(5),
                legend.text = final.subset$Timepoint)
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
#                       VISUALIZE GENE LOCATIONS
#####################################################################################

i <- 2
for (i in 1:length(patients)){
  
  # which patient
  patient <- patients[i]
  print(patient)
  
  # get dataframe subset
  df.subset <- gene.mutation.sample[gene.mutation.sample$Patient == patient, ]
  anno <- c("Timecourse", "EML")
  j <- 1
  
  for (j in 1:length(anno)){
    
    # which anno
    anno.of.interest <- anno[j]
    print(anno.of.interest)
    
    #
    df.anno.subset <- df.subset[df.subset$Annotation == anno.of.interest, ]
    clones <- unique(df.anno.subset$Culture)
    k <- 1
    
    for (k in 1:length(clones)){
      
      # which clone
      clone.of.interest <- clones[k]
      print(clone.of.interest)
      
      #
      final.subset <- df.anno.subset[df.anno.subset$Culture == clone.of.interest, ]
      
      if (anno.of.interest == "Timecourse"){
        
        final.subset$Timepoint <- factor(final.subset$Timepoint, levels = c("101d", "115d", "173d", "190d", "260d", "264d", "404d", "428d", "609d", "705d", "718d", "729d"))
        
      } else if(anno.of.interest == "EML"){
        
        final.subset$Timepoint <- factor(final.subset$Timepoint, levels = c("MID", "LATE"))
        
      }
      
      # reshape it
      matrix <- acast(final.subset, final.subset$Sample ~ final.subset$Type, value.var = "Fraction")
      
      # # order according to mutations
      # final.subset <- final.subset[order(final.subset$total_muts),]
      
      if (anno.of.interest == "EML"){
        
        matrix <- matrix[c(grep("MID", rownames(matrix)),grep("LATE", rownames(matrix))),]
        
        barplot(t(matrix),
                ylab = "Signatures",axes = F,  axisnames=F, ylim = c(0, 1),
                col=c(RColorBrewer::brewer.pal(5, "Set2"),"magenta","firebrick"), cex.lab=1.5, cex.axis=1.5,
                legend.text = colnames(matrix),
                beside = FALSE) # Stacked bars (default)
        
        
        box()
        
      } else {
        
        barplot(t(matrix), ylim = c(0, 1),
                ylab = "Signatures", axes = F,  axisnames=F,
                col=c(RColorBrewer::brewer.pal(5, "Set2"),"magenta","firebrick"), cex.lab=1.5, cex.axis=1.5,
                legend.text = colnames(matrix),
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

# #####################################################################################
# #                       VISUALIZE KATAEGIS
# #####################################################################################

i <- 1
for (i in 1:length(patients)){
  
  # which patient
  patient <- patients[i]
  print(patient)
  
  # get dataframe subset
  df.subset <- kataegis.sample[kataegis.sample$Patient == patient, ]
  anno <- c("Timecourse", "EML")
  j <- 1
  
  for (j in 1:length(anno)){
    
    # which anno
    anno.of.interest <- anno[j]
    print(anno.of.interest)
    
    #
    df.anno.subset <- df.subset[df.subset$Annotation == anno.of.interest, ]
    clones <- unique(df.anno.subset$Culture)
    k <- 1
    
    for (k in 1:length(clones)){
      
      # which clone
      clone.of.interest <- clones[k]
      print(clone.of.interest)
      
      #
      final.subset <- df.anno.subset[df.anno.subset$Culture == clone.of.interest, ]
      
      if (anno.of.interest == "Timecourse"){
        
        final.subset$Timepoint <- factor(final.subset$Timepoint, levels = c("101d", "115d", "173d", "190d", "260d", "264d", "404d", "428d", "609d", "705d", "718d", "729d"))
        
      } else if(anno.of.interest == "EML"){
        
        final.subset$Timepoint <- factor(final.subset$Timepoint, levels = c("MID", "LATE"))
        
      }
      
      if (anno.of.interest == "EML"){
        
        barplot(final.subset$Count, ylim = c(0, 3),
                ylab = "Kataegis", axes = F,  axisnames=F, cex.lab=1.5, cex.axis=1.5, cex.main=1.5,
                col = c("darkblue", "darkmagenta"),
                legend.text = final.subset$Timepoint) # Stacked bars (default)
        box()
        
      } else {
        
        barplot(final.subset$Count,  ylim = c(0, 3), axes = F,
                ylab = "Kataegis", axisnames=F, cex.lab=1.5, cex.axis=1.5, cex.main=1.5,
                col = colfunc(5),
                legend.text = final.subset$Timepoint)
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
  anno <- c("Timecourse", "EML")
  j <- 1
  
  for (j in 1:length(anno)){
    
    # which anno
    anno.of.interest <- anno[j]
    print(anno.of.interest)
    
    #
    df.anno.subset <- df.subset[df.subset$Annotation == anno.of.interest, ]
    clones <- unique(df.anno.subset$Culture)
    k <- 1
    
    for (k in 1:length(clones)){
      
      # which clone
      clone.of.interest <- clones[k]
      print(clone.of.interest)
      
      #
      final.subset <- df.anno.subset[df.anno.subset$Culture == clone.of.interest, ]
      
      if (anno.of.interest == "Timecourse"){
        
        final.subset$Timepoint <- factor(final.subset$Timepoint, levels = c("101d", "115d", "173d", "190d", "260d", "264d", "404d", "428d", "609d", "705d", "718d", "729d"))
        
      }
      
      # reshape it
      matrix <-acast(final.subset, final.subset$Sample ~ final.subset$Type, value.var = "Count")
      
      # # order according to mutations
      # final.subset <- final.subset[order(final.subset$total_muts),]
      
      if (anno.of.interest == "EML"){
        
        matrix <- matrix[c(grep("MID", rownames(matrix)),grep("LATE", rownames(matrix))),]
        
        barplot(t(matrix),
                ylab = "Indel",axes = F,  axisnames=F, ylim = c(0, 800), 
                col = c("darkorange", "purple"), cex.lab=1.5, cex.axis=1.5,
                legend.text = colnames(matrix),
                beside = FALSE) # Stacked bars (default)
        
        
        box()
        
      } else {
        
        barplot(t(matrix), ylim = c(0, 800), 
                ylab = "CNV", axes = F,  axisnames=F,
                col = c("darkorange", "purple"), cex.lab=1.5, cex.axis=1.5,
                legend.text = colnames(matrix),
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
  anno <- c("Timecourse", "EML")
  j <- 1
  
  for (j in 1:length(anno)){
    
    # which anno
    anno.of.interest <- anno[j]
    print(anno.of.interest)
    
    #
    df.anno.subset <- df.subset[df.subset$Annotation == anno.of.interest, ]
    clones <- unique(df.anno.subset$Culture)
    k <- 1
    
    for (k in 1:length(clones)){
      
      # which clone
      clone.of.interest <- clones[k]
      print(clone.of.interest)
      
      final.subset <- df.anno.subset[df.anno.subset$Culture == clone.of.interest, ]
      
      if (anno.of.interest == "Timecourse"){
        
        final.subset$Timepoint <- factor(final.subset$Timepoint, levels = c("101d", "115d", "173d", "190d", "260d", "264d", "404d", "428d", "609d", "705d", "718d", "729d"))
        
      } else if(anno.of.interest == "EML"){
        
        final.subset$Timepoint <- factor(final.subset$Timepoint, levels = c("MID", "LATE"))
        
      }
      
      if (anno.of.interest == "EML"){
        
        barplot(final.subset$Fraction, ylim = c(0, 0.2), main = clone.of.interest,
                ylab = "wGII", axes = F,  axisnames=F, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, width=c(5, 5),
                col = c("darkblue", "darkmagenta"),
                legend.text = final.subset$Timepoint) # Stacked bars (default)
        box()
        
      } else {
        
        barplot(final.subset$Fraction,  ylim = c(0, 0.2), main = clone.of.interest, axes = F,
                ylab = "wGII", axisnames=F, cex.lab=1.5, cex.axis=1.5, cex.main=1.5,
                col =  colfunc(5),
                legend.text = final.subset$Timepoint)
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
#                            VISUALIZE LOH
#####################################################################################

i <- 1
for (i in 1:length(patients)){
  
  # which patient
  patient <- patients[i]
  print(patient)
  
  # get dataframe subset
  df.subset <- loh.data.melt[loh.data.melt$Patient == patient, ]
  anno <- c("Timecourse", "EML")
  j <- 1
  
  for (j in 1:length(anno)){
    
    # which anno
    anno.of.interest <- anno[j]
    print(anno.of.interest)
    
    #
    df.anno.subset <- df.subset[df.subset$Annotation == anno.of.interest, ]
    clones <- unique(df.anno.subset$Culture)
    k <- 1
    
    for (k in 1:length(clones)){
      
      # which clone
      clone.of.interest <- clones[k]
      print(clone.of.interest)
      
      final.subset <- df.anno.subset[df.anno.subset$Culture == clone.of.interest, ]
      
      if (anno.of.interest == "Timecourse"){
        
        final.subset$Timepoint <- factor(final.subset$Timepoint, levels = c("101d", "115d", "173d", "190d", "260d", "264d", "404d", "428d", "609d", "705d", "718d", "729d"))
        
      } else if(anno.of.interest == "EML"){
        
        final.subset$Timepoint <- factor(final.subset$Timepoint, levels = c("MID", "LATE"))
        
      }
      
      if (anno.of.interest == "EML"){
        
        barplot(final.subset$Fraction, ylim = c(0, 0.2), main = clone.of.interest,
                ylab = "LOH", axes = F,  axisnames=F, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, width=c(5, 5),
                col = c("darkblue", "darkmagenta"),
                legend.text = final.subset$Timepoint) # Stacked bars (default)
        box()
        
      } else {
        
        barplot(final.subset$Fraction,  ylim = c(0, 0.2), main = clone.of.interest, axes = F,
                ylab = "LOH", axisnames=F, cex.lab=1.5, cex.axis=1.5, cex.main=1.5,
                col =  colfunc(5),
                legend.text = final.subset$Timepoint)
        box()
        
      }
      
      # add axis only to the first plot
      if (i == 1 & j == 1 & k == 1){
        axis(2, font=1.5)
      }
    }
  }
}

# #####################################################################################
# #                       VISUALIZE COPY NUMBER VARIATION
# #####################################################################################

i <- 1
for (i in 1:length(patients)){
  
  # which patient
  patient <- patients[i]
  print(patient)
  
  # get dataframe subset
  df.subset <- cnv.melt.summary[cnv.melt.summary$Patient == patient, ]
  anno <- c("Timecourse", "EML")
  j <- 2
  
  for (j in 1:length(anno)){
    
    # which anno
    anno.of.interest <- anno[j]
    print(anno.of.interest)
    
    #
    df.anno.subset <- df.subset[df.subset$Annotation == anno.of.interest, ]
    clones <- unique(df.anno.subset$Culture)
    k <- 1
    
    for (k in 1:length(clones)){
      
      # which clone
      clone.of.interest <- clones[k]
      print(clone.of.interest)
      
      #
      final.subset <- df.anno.subset[df.anno.subset$Culture == clone.of.interest, ]
      
      if (anno.of.interest == "Timecourse"){
        
        final.subset$Timepoint <- factor(final.subset$Timepoint, levels = c("101d", "115d", "173d", "190d", "260d", "264d", "404d", "428d", "609d", "705d", "718d", "729d"))
        
      }
      
      # reshape it
      matrix <-acast(final.subset, final.subset$Sample ~ final.subset$Type, value.var = "Count")
      
      # # order according to mutations
      # final.subset <- final.subset[order(final.subset$total_muts),]
      
      if (anno.of.interest == "EML"){
        
        matrix <- matrix[c(grep("MID", rownames(matrix)),grep("LATE", rownames(matrix))),]
        
        barplot(t(matrix),
                ylab = "CNV",axes = F,  axisnames=F, ylim = c(0, 10),
                col = c("red", "darkblue"), cex.lab=1.5, cex.axis=1.5,
                legend.text = colnames(matrix),
                beside = FALSE) # Stacked bars (default)
        
        
        box()
        
      } else {
        
        barplot(t(matrix), ylim = c(0, 10),
                ylab = "CNV", axes = F,  axisnames=F,
                col = c("red", "darkblue"), cex.lab=1.5, cex.axis=1.5,
                legend.text = colnames(matrix),
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
#                       VISUALIZE STRUCTURAL VARIANTS
#####################################################################################

i <- 1
for (i in 1:length(patients)){
  
  # which patient
  patient <- patients[i]
  print(patient)
  
  # get dataframe subset
  df.subset <- sv.sample[sv.sample$Patient == patient, ]
  anno <- c("Timecourse", "EML")
  j <- 1
  
  for (j in 1:length(anno)){
    
    # which anno
    anno.of.interest <- anno[j]
    print(anno.of.interest)
    
    #
    df.anno.subset <- df.subset[df.subset$Annotation == anno.of.interest, ]
    clones <- unique(df.anno.subset$Culture)
    k <- 1
    
    for (k in 1:length(clones)){
      
      # which clone
      clone.of.interest <- clones[k]
      print(clone.of.interest)
      
      #
      final.subset <- df.anno.subset[df.anno.subset$Culture == clone.of.interest, ]
      
      if (anno.of.interest == "Timecourse"){
        
        final.subset$Timepoint <- factor(final.subset$Timepoint, levels = c("101d", "115d", "173d", "190d", "260d", "264d", "404d", "428d", "609d", "705d", "718d", "729d"))
        
      }
      
      # reshape it
      matrix <-acast(final.subset, final.subset$Sample ~ final.subset$Type, value.var = "Count")
      matrix[is.na(matrix)] <- "0"
      
      # # order according to mutations
      # final.subset <- final.subset[order(final.subset$total_muts),]
      
      if (anno.of.interest == "EML"){
        
        matrix <- matrix[c(grep("MID", rownames(matrix)),grep("LATE", rownames(matrix))),]
        
        barplot(t(matrix),
                ylab = "SV",axes = F,  axisnames=F, ylim = c(0, 100), 
                col = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"), cex.lab=1.5, cex.axis=1.5,
                legend.text = colnames(matrix),
                beside = FALSE) # Stacked bars (default)
        
        
        box()
        
      } else {
        
        barplot(t(matrix), ylim = c(0, 100), 
                ylab = "SV", axes = F,  axisnames=F,
                col = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"), cex.lab=1.5, cex.axis=1.5,
                legend.text = colnames(matrix),
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
#                       VISUALIZE PROPORTION OF SVS OVER TIME
#####################################################################################

i <- 1
for (i in 1:length(patients)){
  
  # which patient
  patient <- patients[i]
  print(patient)
  
  # get dataframe subset
  df.subset <- sv.prop[sv.prop$Patient == patient, ]
  anno <- c("Timecourse", "EML")
  j <- 1
  
  for (j in 1:length(anno)){
    
    # which anno
    anno.of.interest <- anno[j]
    print(anno.of.interest)
    
    #
    df.anno.subset <- df.subset[df.subset$Annotation == anno.of.interest, ]
    clones <- unique(df.anno.subset$Culture)
    k <- 1
    
    for (k in 1:length(clones)){
      
      # which clone
      clone.of.interest <- clones[k]
      print(clone.of.interest)
      
      #
      final.subset <- df.anno.subset[df.anno.subset$Culture == clone.of.interest, ]
      
      if (anno.of.interest == "Timecourse"){
        
        final.subset$Timepoint <- factor(final.subset$Timepoint, levels = c("101d", "115d", "173d", "190d", "260d", "264d", "404d", "428d", "609d", "705d", "718d", "729d"))
        
      } else if(anno.of.interest == "EML"){
        
        final.subset$Timepoint <- factor(final.subset$Timepoint, levels = c("MID", "LATE"))
        
      }
      
      # reshape it
      matrix <- acast(final.subset, final.subset$Sample ~ final.subset$Type, value.var = "Fraction")
      
      # # order according to mutations
      # final.subset <- final.subset[order(final.subset$total_muts),]
      
      if (anno.of.interest == "EML"){
        
        matrix <- matrix[c(grep("MID", rownames(matrix)),grep("LATE", rownames(matrix))),]
        
        barplot(t(matrix),
                ylab = "SV frequency",axes = F,  axisnames=F, ylim = c(0, 1),
                col = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"), cex.lab=1.5, cex.axis=1.5,
                legend.text = colnames(matrix),
                beside = FALSE) # Stacked bars (default)
        
        
        box()
        
      } else {
        
        barplot(t(matrix), ylim = c(0, 1),
                ylab = "SV frequency", axes = F,  axisnames=F,
                col = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"), cex.lab=1.5, cex.axis=1.5,
                legend.text = colnames(matrix),
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

