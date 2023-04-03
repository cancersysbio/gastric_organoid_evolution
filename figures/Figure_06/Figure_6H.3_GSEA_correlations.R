#############################################################################################################################
##                                                                                                                      
##  SCRIPT TO CREATE SPEARMAN CORRELATION FOR GSEA OUTPUT
##                                                                                                                      
##  
##  Author: Kasper Karlsson
##
##                                                                                                                      
############################################################################################################################

library(tidyr)
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)


setwd(file.path("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_6/Figure_6H"))

D2C1R1 <- read.table("D2C1R1/GSEA_compiled_all_direction_D2C1R1.txt",sep="\t",header=TRUE)
D2C1R2 <- read.table("D2C1R2/GSEA_compiled_all_direction_D2C1R2.txt",sep="\t",header=TRUE)
D2C2R2 <- read.table("D2C2R2/GSEA_compiled_all_direction_D2C2R2.txt",sep="\t",header=TRUE)
D3C2R1 <- read.table("D3C2R1/GSEA_compiled_all_direction_D3C2R1.txt",sep="\t",header=TRUE)
EML <- read.table("EML/GSEA_compiled_all_direction_EML.txt",sep="\t",header=TRUE)

D2C1R1$Timepoint <- as.character(D2C1R1$Timepoint)
D2C1R2$Timepoint <- as.character(D2C1R2$Timepoint)
D2C2R2$Timepoint <- as.character(D2C2R2$Timepoint)
D3C2R1$Timepoint <- as.character(D3C2R1$Timepoint)
EML$Culture <- paste0(EML$Culture,"_",EML$Timepoint)
combined <- bind_rows(D2C1R1,D2C1R2,D2C2R2,D3C2R1,EML, .id = "column_label")

combined2 <- subset(combined, Timepoint != "315") # Remove latest time point for ECB replicates for simplification
all_small <- select(combined2, Pathway, Culture, Observed_Score)

Pathways <- c("HALLMARK_MYC_TARGETS_V1","HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_E2F_TARGETS","HALLMARK_HYPOXIA","HALLMARK_OXIDATIVE_PHOSPHORYLATION",
              "HALLMARK_G2M_CHECKPOINT","HALLMARK_APOPTOSIS","HALLMARK_MTORC1_SIGNALING","HALLMARK_P53_PATHWAY","HALLMARK_MYC_TARGETS_V2")

all_small2 <- subset(all_small, Pathway %in% Pathways)

samples <- c("D1C1_LATE","D1C2_LATE","D1C3_LATE","D2C2_LATE","D2C3_LATE","D3C2_LATE",  
             "D2C1R1_0a","D2C1R2_1a","D2C2R2_0a","D3C2R1_0a",                           
             "D2C1R1_0b","D2C1R1_1b","D2C1R1_2","D2C1R1_4",                          
             "D2C1R2_0a","D2C1R2_0b","D2C1R2_1b","D2C1R2_2","D2C1R2_4",                
             "D2C2R2_0b","D2C2R2_1","D2C2R2_11","D2C2R2_13","D2C2R2_39","D2C2R2_56","D2C2R2_6","D2C2R2_9",
             "D3C2R1_0b","D3C2R1_1","D3C2R1_2a","D3C2R1_2b","D3C2R1_2c","D3C2R1_3","D3C2R1_4","D3C2R1_7","D3C2R1_8")

all_small3 <- subset(all_small2, Culture %in% samples)
all_small_long <- reshape(all_small3, idvar = "Pathway", timevar = "Culture", direction = "wide")
rownames(all_small_long) <- all_small_long$Pathway
all_small_long$Pathway <- NULL
cc = cor(all_small_long, method = "spearman")

pdf("Complex_heatmap_full.pdf", width = 12, height = 12)

Heatmap(cor(all_small_long,method="spearman"), name = "cor", 
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), 
        #right_annotation = ha_row, 
        show_row_names = TRUE, show_column_names = FALSE, row_dend_side = NULL, 
        show_column_dend = TRUE, column_title = "pairwise correlation between samples",
        heatmap_legend_param = list(title = "Correlation"))
dev.off()






