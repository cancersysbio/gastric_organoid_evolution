#############################################################################################################################
##                                                                                                                      
##  GET DIFFERENTIAL GENE EXPRESSION FROM SEURAT OBJECT FOR DONOR 1
##                                                                                                                      
##  
##  Author: Kasper Karlsson
##
##                                                                                                                      
############################################################################################################################


### USE R version 4.0

library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(tidyr)


### Load Seurat object
setwd(file.path("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Data_objects/EML_Seq_D1"))

seq_EML_D1 <- readRDS("EML_Seq_D1_seurat_obj.rds")
head(seq_EML_D1@meta.data)
Idents(seq_EML_D1) <- seq_EML_D1@meta.data$orig.ident
table(seq_EML_D1@meta.data$orig.ident)


########## DIFFERENTIAL GENE EXPRESSION


### D1
markers_D1_C1_Early_vs_WT <- FindMarkers(seq_EML_D1, ident.1="D1_C1_Early", ident.2="D1_WT", logfc.threshold = 0)
markers_D1_C2_Early_vs_WT <- FindMarkers(seq_EML_D1, ident.1="D1_C2_Early", ident.2="D1_WT", logfc.threshold = 0)
markers_D1_C3_Early_vs_WT <- FindMarkers(seq_EML_D1, ident.1="D1_C3_Early", ident.2="D1_WT", logfc.threshold = 0)

markers_D1_C1_Mid_vs_WT <- FindMarkers(seq_EML_D1, ident.1="D1_C1_Mid", ident.2="D1_WT", logfc.threshold = 0)
markers_D1_C2_Mid_vs_WT <- FindMarkers(seq_EML_D1, ident.1="D1_C2_Mid", ident.2="D1_WT", logfc.threshold = 0)
markers_D1_C3_Mid_vs_WT <- FindMarkers(seq_EML_D1, ident.1="D1_C3_Mid", ident.2="D1_WT", logfc.threshold = 0)

markers_D1_C1_Late_vs_WT <- FindMarkers(seq_EML_D1, ident.1="D1_C1_Late", ident.2="D1_WT", logfc.threshold = 0)
markers_D1_C2_Late_vs_WT <- FindMarkers(seq_EML_D1, ident.1="D1_C2_Late", ident.2="D1_WT", logfc.threshold = 0)
markers_D1_C3_Late_vs_WT <- FindMarkers(seq_EML_D1, ident.1="D1_C3_Late", ident.2="D1_WT", logfc.threshold = 0)

write.table(markers_D1_C1_Early_vs_WT,quote=FALSE,paste0("DEGs_all/all_markers_D1_C1_Early_vs_WT.txt"), sep="\t")
write.table(markers_D1_C2_Early_vs_WT,quote=FALSE,paste0("DEGs_all/all_markers_D1_C2_Early_vs_WT.txt"), sep="\t")
write.table(markers_D1_C3_Early_vs_WT,quote=FALSE,paste0("DEGs_all/all_markers_D1_C3_Early_vs_WT.txt"), sep="\t")
write.table(markers_D1_C1_Mid_vs_WT,quote=FALSE,paste0("DEGs_all/all_markers_D1_C1_Mid_vs_WT.txt"), sep="\t")
write.table(markers_D1_C2_Mid_vs_WT,quote=FALSE,paste0("DEGs_all/all_markers_D1_C2_Mid_vs_WT.txt"), sep="\t")
write.table(markers_D1_C3_Mid_vs_WT,quote=FALSE,paste0("DEGs_all/all_markers_D1_C3_Mid_vs_WT.txt"), sep="\t")
write.table(markers_D1_C1_Late_vs_WT,quote=FALSE,paste0("DEGs_all/all_markers_D1_C1_Late_vs_WT.txt"), sep="\t")
write.table(markers_D1_C2_Late_vs_WT,quote=FALSE,paste0("DEGs_all/all_markers_D1_C2_Late_vs_WT.txt"), sep="\t")
write.table(markers_D1_C3_Late_vs_WT,quote=FALSE,paste0("DEGs_all/all_markers_D1_C3_Late_vs_WT.txt"), sep="\t")

markers_D1_C1_Mid_vs_Early <- FindMarkers(seq_EML_D1, ident.1="D1_C1_Mid", ident.2="D1_C1_Early", logfc.threshold = 0)
markers_D1_C2_Mid_vs_Early <- FindMarkers(seq_EML_D1, ident.1="D1_C2_Mid", ident.2="D1_C2_Early", logfc.threshold = 0)
markers_D1_C3_Mid_vs_Early <- FindMarkers(seq_EML_D1, ident.1="D1_C3_Mid", ident.2="D1_C3_Early", logfc.threshold = 0)

markers_D1_C1_Late_vs_Early <- FindMarkers(seq_EML_D1, ident.1="D1_C1_Late", ident.2="D1_C1_Early", logfc.threshold = 0)
markers_D1_C2_Late_vs_Early <- FindMarkers(seq_EML_D1, ident.1="D1_C2_Late", ident.2="D1_C2_Early", logfc.threshold = 0)
markers_D1_C3_Late_vs_Early <- FindMarkers(seq_EML_D1, ident.1="D1_C3_Late", ident.2="D1_C3_Early", logfc.threshold = 0)


write.table(markers_D1_C1_Mid_vs_Early,quote=FALSE,paste0("DEGs_all/vsEarly/all_markers_D1_C1_Mid_vs_Early.txt"), sep="\t")
write.table(markers_D1_C2_Mid_vs_Early,quote=FALSE,paste0("DEGs_all/vsEarly/all_markers_D1_C2_Mid_vs_Early.txt"), sep="\t")
write.table(markers_D1_C3_Mid_vs_Early,quote=FALSE,paste0("DEGs_all/vsEarly/all_markers_D1_C3_Mid_vs_Early.txt"), sep="\t")
write.table(markers_D1_C1_Late_vs_Early,quote=FALSE,paste0("DEGs_all/vsEarly/all_markers_D1_C1_Late_vs_Early.txt"), sep="\t")
write.table(markers_D1_C2_Late_vs_Early,quote=FALSE,paste0("DEGs_all/vsEarly/all_markers_D1_C2_Late_vs_Early.txt"), sep="\t")
write.table(markers_D1_C3_Late_vs_Early,quote=FALSE,paste0("DEGs_all/vsEarly/all_markers_D1_C3_Late_vs_Early.txt"), sep="\t")


markers_D1_C1_Late_vs_Mid <- FindMarkers(seq_EML_D1, ident.1="D1_C1_Late", ident.2="D1_C1_Mid", logfc.threshold = 0)
markers_D1_C2_Late_vs_Mid <- FindMarkers(seq_EML_D1, ident.1="D1_C2_Late", ident.2="D1_C2_Mid", logfc.threshold = 0)
markers_D1_C3_Late_vs_Mid <- FindMarkers(seq_EML_D1, ident.1="D1_C3_Late", ident.2="D1_C3_Mid", logfc.threshold = 0)

write.table(markers_D1_C1_Late_vs_Mid,quote=FALSE,paste0("DEGs_all/vsEML/all_markers_D1_C1_Late_vs_Mid.txt"), sep="\t")
write.table(markers_D1_C2_Late_vs_Mid,quote=FALSE,paste0("DEGs_all/vsEML/all_markers_D1_C2_Late_vs_Mid.txt"), sep="\t")
write.table(markers_D1_C3_Late_vs_Mid,quote=FALSE,paste0("DEGs_all/vsEML/all_markers_D1_C3_Late_vs_Mid.txt"), sep="\t")




