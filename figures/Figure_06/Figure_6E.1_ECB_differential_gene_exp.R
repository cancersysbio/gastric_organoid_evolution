#############################################################################################################################
##                                                                                                                      
##  SCRIPT TO GET ECB SPECIFIC DIFFERENTIAL GENE EXPRESSION
##                                                                                                                      
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
seurat_seq8 <- readRDS("Sequencing8_D2C2T2R2_seurat_obj.rds")



######### UPDATE CELL ANNOTATION WITH SUBCLONE 0A

Idents(seurat_seq8) <- seurat_seq8@meta.data$ECB_RG
RG_0a <- scan("RG0_4q-20q+.txt", character(), quote = "") ### List of cells that are 4q- and 20q+ (from inferCNV)
seurat_seq8 <- SetIdent(seurat_seq8, cells = RG_0a, value = "0a")

seurat_seq8@meta.data$ECB_updated <- Idents(seurat_seq8)
unique(seurat_seq8@meta.data$ECB_updated)



### Differential Gene Expression based on RG


rgs <- c("0a","0","1","39","13","56","11","6","9")   



for (r in rgs) {
  print (r)
  markers <- FindMarkers(seurat_seq8, ident.1=r,  logfc.threshold = 0)
  write.table(markers,quote=FALSE,paste0("DEGs_allGenes/markers_RG_",r,".txt"), sep="\t")
}

markers <- FindMarkers(seurat_seq8, ident.1="0a", ident.2 = "0", logfc.threshold = 0)
write.table(markers,quote=FALSE,paste0("DEGs_allGenes_againstRG0/markers_RG_0a_vs_0b.txt"), sep="\t")


