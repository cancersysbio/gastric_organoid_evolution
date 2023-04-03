#############################################################################################################################
##                                                                                                                      
##  GET DIFFERENTIAL GENE EXPRESSION FROM SEURAT OBJECT FOR DONOR 2
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


version

### Load Seurat object
setwd(file.path("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Data_objects/EML_Seq_D2"))

seq_EML_D2 <- readRDS("EML_Seq_D2_seurat_obj.rds")

head(seq_EML_D2@meta.data)
Idents(seq_EML_D2) <- seq_EML_D2@meta.data$orig.ident
table(seq_EML_D2@meta.data$orig.ident)


########## DIFFERENTIAL GENE EXPRESSION

### D2

markers_D2_C2_Early_vs_WT <- FindMarkers(seq_EML_D2, ident.1="D2_C2_Early", ident.2="D2_WT", logfc.threshold = 0)
markers_D2_C3_Early_vs_WT <- FindMarkers(seq_EML_D2, ident.1="D2_C3_Early", ident.2="D2_WT", logfc.threshold = 0)

markers_D2_C1_Mid_vs_WT <- FindMarkers(seq_EML_D2, ident.1="D2_C1_Mid", ident.2="D2_WT", logfc.threshold = 0)
markers_D2_C2_Mid_vs_WT <- FindMarkers(seq_EML_D2, ident.1="D2_C2_Mid", ident.2="D2_WT", logfc.threshold = 0)
markers_D2_C3_Mid_vs_WT <- FindMarkers(seq_EML_D2, ident.1="D2_C3_Mid", ident.2="D2_WT", logfc.threshold = 0)

markers_D2_C2_Late_vs_WT <- FindMarkers(seq_EML_D2, ident.1="D2_C2_Late", ident.2="D2_WT", logfc.threshold = 0)
markers_D2_C3_Late_vs_WT <- FindMarkers(seq_EML_D2, ident.1="D2_C3_Late", ident.2="D2_WT", logfc.threshold = 0)

write.table(markers_D2_C2_Early_vs_WT,quote=FALSE,paste0("DEGs_all/all_markers_D2_C2_Early_vs_WT.txt"), sep="\t")
write.table(markers_D2_C3_Early_vs_WT,quote=FALSE,paste0("DEGs_all/all_markers_D2_C3_Early_vs_WT.txt"), sep="\t")
write.table(markers_D2_C1_Mid_vs_WT,quote=FALSE,paste0("DEGs_all/all_markers_D2_C1_Mid_vs_WT.txt"), sep="\t")
write.table(markers_D2_C2_Mid_vs_WT,quote=FALSE,paste0("DEGs_all/all_markers_D2_C2_Mid_vs_WT.txt"), sep="\t")
write.table(markers_D2_C3_Mid_vs_WT,quote=FALSE,paste0("DEGs_all/all_markers_D2_C3_Mid_vs_WT.txt"), sep="\t")
write.table(markers_D2_C2_Late_vs_WT,quote=FALSE,paste0("DEGs_all/all_markers_D2_C2_Late_vs_WT.txt"), sep="\t")
write.table(markers_D2_C3_Late_vs_WT,quote=FALSE,paste0("DEGs_all/all_markers_D2_C3_Late_vs_WT.txt"), sep="\t")

markers_D2_C2_Mid_vs_Early <- FindMarkers(seq_EML_D2, ident.1="D2_C2_Mid", ident.2="D2_C2_Early", logfc.threshold = 0)
markers_D2_C2_Late_vs_Early <- FindMarkers(seq_EML_D2, ident.1="D2_C2_Late", ident.2="D2_C2_Early", logfc.threshold = 0)
markers_D2_C3_Mid_vs_Early <- FindMarkers(seq_EML_D2, ident.1="D2_C3_Mid", ident.2="D2_C3_Early", logfc.threshold = 0)
markers_D2_C3_Late_vs_Early <- FindMarkers(seq_EML_D2, ident.1="D2_C3_Late", ident.2="D2_C3_Early", logfc.threshold = 0)
markers_D2_C2_Late_vs_Mid <- FindMarkers(seq_EML_D2, ident.1="D2_C2_Late", ident.2="D2_C2_Mid", logfc.threshold = 0)
markers_D2_C3_Late_vs_Mid <- FindMarkers(seq_EML_D2, ident.1="D2_C3_Late", ident.2="D2_C3_Mid", logfc.threshold = 0)

write.table(markers_D2_C2_Mid_vs_Early,quote=FALSE,paste0("DEGs_all/vsEarly/all_markers_D2_C2_Mid_vs_Early.txt"), sep="\t")
write.table(markers_D2_C2_Late_vs_Early,quote=FALSE,paste0("DEGs_all/vsEarly/all_markers_D2_C2_Late_vs_Early.txt"), sep="\t")
write.table(markers_D2_C3_Mid_vs_Early,quote=FALSE,paste0("DEGs_all/vsEarly/all_markers_D2_C3_Mid_vs_Early.txt"), sep="\t")
write.table(markers_D2_C3_Late_vs_Early,quote=FALSE,paste0("DEGs_all/vsEarly/all_markers_D2_C3_Late_vs_Early.txt"), sep="\t")
write.table(markers_D2_C2_Late_vs_Mid,quote=FALSE,paste0("DEGs_all/vsEML/all_markers_D2_C2_Late_vs_Mid.txt"), sep="\t")
write.table(markers_D2_C3_Late_vs_Mid,quote=FALSE,paste0("DEGs_all/vsEML/all_markers_D2_C3_Late_vs_Mid.txt"), sep="\t")



