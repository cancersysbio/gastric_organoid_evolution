#############################################################################################################################
##                                                                                                                      
##  SINGLE CELL RNA SEQUENCING ANALYSIS
##  FOR FIGURE 6
##                                                                                                                      
##  Date: 27 DECEMBER 2021                                                                                                                   
##  
##  Author: Kasper Karlsson
##
##                                                                                                                      
############################################################################################################################

library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(gghighlight)


### Load Seurat object
setwd(file.path("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Data_objects/Seq8_D2C2R2T2"))
seurat_D2C2R2 <- readRDS("Sequencing8_D2C2T2R2_seurat_obj.rds")

### Check number of cells and number of cells per barcode
Idents(seurat_D2C2R2) <- seurat_D2C2R2@meta.data$ECB_RG
seurat_D2C2R2
write.table(sort(table(seurat_D2C2R2@meta.data$ECB_RG),decreasing = TRUE),file="scRNA_cells_per_barcode.csv",sep=",")

### Update Seurat object with inferCNV results
RG_0a <- scan("RG0_4q-20q+.txt", character(), quote = "") ### List of cells that are 4q- and 20q+ (from inferCNV)
seurat_D2C2R2 <- SetIdent(seurat_D2C2R2, cells = RG_0a, value = "0a")
seurat_D2C2R2@meta.data$ECB_updated <- Idents(seurat_D2C2R2)

### DIFFERENTIAL GENE EXPRESSION BASED ON RG (FOR FIGURES 6D-6F)
rgs <- c("0a","0","1","39","13","56","11","6","9")   


### EXTRACT TOP DIFFERENTIALLY EXPRESSED GENES IN STAD FROM GEPIA (FOR FIGURE 6E)
gepia <- read.csv("table_degenes_STAD.txt",header=TRUE,sep="\t")
gepia <- gepia[order(-gepia$Log2.Fold.Change.),] 
gepiaT20.genes <- as.vector(gepia$Gene.Symbol[1:20])


### HEATMAP OF TOP GASTRIC CANCER EXPRESSED GENES (ACCORDING TO TCGA - THROUGH GEPIA FOR STAD)
subset_topBC <- subset(seurat_D2C2R2,idents = c("0a","0","1","39","13","56","11","6","9")) ### 11 cells
seurat.obj.big <- subset_topBC
levels(seurat.obj.big) <- c("0a","0","1","39","13","56","11","6","9")

DotPlot(seurat.obj.big, features = gepiaT20.genes, cols = c("lightgrey", "darkred"), scale = TRUE, col.min = 0, col.max = 3) +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=10, face = "bold")) + 
  theme(legend.text = element_text(colour="black", size=10, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=12),
        axis.text.x = element_text(face="bold", color="black", size=12, angle = 45, hjust = 1))


ggsave(filename = "DotPlot_specific_genes/D2C2R2_top_genpia_genes.pdf", width = 4, height = 6, units = "in", device='pdf')
