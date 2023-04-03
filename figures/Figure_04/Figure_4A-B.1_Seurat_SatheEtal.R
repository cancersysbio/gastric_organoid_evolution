#############################################################################################################################
##                                                                                                                      
##  Analysis of single-cell RNA-seq data Sathe et al., 2020, Clinical Cancer Research using Seurat
##                                                                                                                      
##  Date: 23 July 2020                                                                                                                    
##  
##  Author: Moritz Przybilla                                                                                                                    
##                                                                                                                      
#############################################################################################################################

## SEURAT
#' Seurat is an analysis tool developed by Rahul Satija's Lab at the New York Genome Center in Lower Mannhattan. This script is
#' based on their featured Seurat introductions which can be found [here](https://satijalab.org/seurat/). There are various vignettes
#' for distinct analytical workflows and procedures provided. The following code is a breakdown of what you can find there. 
#' In case of questions and unclear points in this Markdown document, it might be best to look into their original vignettes for further descriptions. 
#' The `Seurat` version used in this analysis is `Seurat v3.0` upwards. 
#' Also, it might be beneficial to check out their publication [Stuart, Butler et al., 2019, Cell](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8).
#' The first part of this markdown file will focus on single sample analysis. 
#' The fourth and last part will then explore differential gene expression analysis and further visualization of the data. 

#############################################################################################################################

# clear workspace
rm(list = ls())
set.seed(14) # set the seed of random number generator 

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("tidyverse", "patchwork", "Seurat", "Matrix", "biomaRt", "scater", "DoubletFinder", "viridis")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages, repos = "http://cran.us.r-project.org")

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

# set up functions which are used 
`%notin%` <- Negate(`%in%`)

#############################################################################################################################

## READ IN YOUR DATA OF INTEREST
#' Here, we read in all the 10x Genomics matrices, that belong to the set of samples we want to merge. 
#' Merging data is only recommended if you expect lower batch effects than in for instance the analysis of datasets from distinct
#' laboratories and technologies. Here, we assume that we have minimal batch effects in our data from the same set of experiments
#' and technologies. For the dataset from Zhang et al., 2019, Cell we use the raw count matrices from pre-processed Seurat objects. 

## SATHE ET AL
c.dir <- "/labs/ccurtis2/mjprzy/scRNA_analysis/SatheEtal2020_scRNA_gastric/filtered_data"

# define output directory
o.dir <- "/labs/ccurtis2/mjprzy/scRNA_analysis/hashECB_data_freeze"

# define a sample id for the merged object
sample.tmp <- "merged_data_SatheEtal"

# create sample-specific output directory with a "plots" subfolder
dir.create(paste0(o.dir, "/" , sample.tmp))
dir.create(paste0(o.dir, "/" , sample.tmp, "/plots"))

############################################################################
##               Read in the Seurat objects for SATHE
############################################################################

## SATHE ET AL
sample.list <- list.files(c.dir, full.names = T, include.dirs = T)
sample.list <- sample.list[-grep("pbmc", sample.list)]

# get sample ids for ZHANG AND SATHE
sample.ids <- str_split_fixed(sample.list, "/", 10)[,8]

## SATHE ET AL - 10X MATRICES
sample.list <- lapply(sample.list, Read10X)

# list to store seurat.objects in
seurat.obj.list <- list()

# iterate over seurat objects and extract raw count matrices
for (i in 1:length(sample.list)){
  
  # extract count matrix
  c.mtx <- sample.list[[i]]
  
  # create new seurat object
  seurat.obj <- CreateSeuratObject(counts = c.mtx, 
                                   project = paste0(sample.ids[i], "_GC"), 
                                   min.cells = 3, 
                                   min.features = 200)
  
  seurat.obj.list[[i]] <- seurat.obj
  
}

# combined the seurat objects together
seurat.obj.big <- merge(seurat.obj.list[[1]], y = seurat.obj.list[c(2:length(seurat.obj.list))],
                        add.cell.ids = sample.ids, project = "SatheEtal2020_scRNA_gastric", merge.data = T)

# save identity
seurat.obj.big$sample_ident <- seurat.obj.big$orig.ident

# make new metadata column with sample type
sample.type <- str_split_fixed(seurat.obj.big@meta.data$orig.ident, "_GC", 2)[,1]
sample.type[grep("n1", sample.type)] <- "Normal"
sample.type[grep("n2", sample.type)] <- "Normal"
sample.type[grep("t1", sample.type)] <- "Tumor"
sample.type[grep("t2", sample.type)] <- "Tumor"
seurat.obj.big@meta.data$sample_type <- sample.type

# check colnames
head(colnames(seurat.obj.big))
tail(colnames(seurat.obj.big))

# check amount of cells
table(seurat.obj.big$orig.ident)

############################################################################
##                     Check out the raw matrix first
############################################################################

c.mtx <- seurat.obj.big@assays$RNA@counts

# calculate counts per cell
counts.per.cell <- Matrix::colSums(c.mtx)

# calculate nubmer of genes per cell
counts.per.gene <- Matrix::rowSums(c.mtx)

# caluclate genes per cell
genes.per.cell <- Matrix::colSums(c.mtx>0) # count gene only if it has non-zero reads mapped.

# calculate cells per genes
cells.per.gene <- Matrix::rowSums(c.mtx>0) # only count cells where the gene is expressed

# plot histograms looking at the calculated metrics
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_preSeurat_QCplots.pdf"), width = 20, height = 20)
hist(log10(counts.per.cell+1),main='counts per cell',col='wheat')
hist(log10(genes.per.cell+1), main='genes per cell', col='wheat')
plot(counts.per.cell, genes.per.cell, log='xy', col='wheat')
title('counts vs genes per cell')
hist(log10(counts.per.gene+1), main='counts per gene', col='wheat')
dev.off()

# plot each cell ranked by their number of genes detected per cell
# represents distribution of library complexity in the sequencing run
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_library_complexity.pdf"), width = 20, height = 20)
plot(sort(genes.per.cell), xlab='cell', log='y', main='genes per cell (ordered)')
dev.off()

#############################################################################################################################

## DATA PRE-PROCESSING
#' Now, we want to create a `Seurat object`, which we can further use to perform quality control and biological meaningful analysis.
#' Upon creation of the `Seurat object`, the number of genes and UMIs (nGene and nUMI) is calculated automatically for every object.
#' For non-UMI data, nUMI represents the sum of the non-normalized values within a cell. In the next steps, we further calculate the
#' percentage of `mitochondrial and ribosomal genes`, as well as the number of `housekeeping genes` per cell. This is adapted from an
#' older Seurat version, where this was done manually. In particular, these HK genes reflect commomn processes active in a cell and
#' hence provide a good global quality measure. The list of HK genes is adapted from [Tirosh et al., 2016, Nature](https://pubmed.ncbi.nlm.nih.gov/27806376/).

message("Start pre-processing..")

# Calculate the percentage of mitochondrial genes
seurat.obj.big[["percent.mt"]] <- PercentageFeatureSet(seurat.obj.big, pattern = "^MT-")

# Calculate percent ribosomal genes
seurat.obj.big[["percent.ribo"]] <- PercentageFeatureSet(seurat.obj.big, pattern = "^RP")

# the quality metrices in the seurat object are stored in the metadata
head(seurat.obj.big@meta.data, 5)

# make a FeatureScatter plot to visualize feature-feature relationships (can be used generally for all columns in metadata)
plot1 <- FeatureScatter(seurat.obj.big, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat.obj.big, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# save plot 1 and plot2 together
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_FeatureScatterplots.pdf"), width = 20, height = 20)
print(plot1 + plot2)
dev.off()

# Load the the list of house keeping genes
hkgenes <- read.table("/labs/ccurtis2/mjprzy/scRNA_analysis/housekeeping_genes_tirosh.txt")
hkgenes <- as.vector(hkgenes$V1)

# remove hkgenes that were not found
hkgenes.found <- which(toupper(rownames(seurat.obj.big@assays$RNA@counts)) %in% hkgenes)

# calculate the number of housekeeping genes per cell and add it as metadata
n.expressed.hkgenes <- Matrix::colSums(seurat.obj.big@assays$RNA@counts[hkgenes.found, ] > 0)
seurat.obj.big <- AddMetaData(object = seurat.obj.big, metadata = n.expressed.hkgenes, col.name = "n.exp.hkgenes")

# plot additional QC
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_QC_plots.pdf"), width = 20, height = 20)
print(VlnPlot(object = seurat.obj.big, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "n.exp.hkgenes"), ncol = 5))
dev.off()

#############################################################################################################################

## REMOVE CELLS NOT MATCHING THE QC THRESHOLDS
#' Once this step is implemented, we should visually examine the quality of the dataset. Following this, we will
#' remove cells based on the calculated QC metrics. However, it is recommended to use very conservative thresholds,
#' as we do not want to remove valuable biological insights. 

# subset the seurat object to the QC-passing cells
seurat.obj.big <- subset(seurat.obj.big, subset = nFeature_RNA > 500 & percent.mt < 20 & n.exp.hkgenes > 55 & percent.ribo < 40)

# plot post QC
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_postQC_plots.pdf"), width = 20, height = 20)
print(VlnPlot(object = seurat.obj.big, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "n.exp.hkgenes"), ncol = 5))
dev.off()

# check amount of cells after QC
table(seurat.obj.big$orig.ident)

#############################################################################################################################

## NORMALIZATION WITH SCTRANSFORM
#' Now that we are having our quality controlled cells in the dataset, we can implement the next step - normalization of the data.
#' The Seurat workflow previously featured a log-normalization, now recommending to use an alternative method, called `sctransform`.
#' In contrast to the previously implemented normalization method, `sctransform` uses a regularied negative binomial regression
#' for the normalization and variance stabilization. It does so by modeling the expression of each gene individually, then grouping
#' them together based on similarity in terms of the determined model. The regression model can incorporate various independent
#' variables, i.e. the sequencing depth, mitochondrial genes or the difference in cell cycle marker expression. For detailed information
#' have a look at the [Github repository](https://github.com/ChristophH/sctransform) or the
#' [publication](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1).
#' 
#' In this step, we will again use the list of cell cycle genes from [Tirosh et al., 2015, Nature](), which is loaded automatically
#' within the Seurat package. We can segregate this list into markers of **G2/M phase** and markers of **S phase**. Based on this we
#' can determine the signals separating non-cycling cells and cycling cells. Using this knowledge, we can provide the cell cycle difference 
#' as an independent variable to `sctransform`, where we will regress out the differences in cell cycle phase amongst proliferating cells 
#' (which are often uninteresting).

# get cell cycle marker genes from Tirosh et al.
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# normalize data with SCTransform() before applying the CellCycleScoring
seurat.obj.big <- SCTransform(seurat.obj.big, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA'))

# perform cell cycle analysis (make sure to specify the "assay" parameter)
seurat.obj.big <- CellCycleScoring(seurat.obj.big, s.features = s.genes, g2m.features = g2m.genes, assay = 'SCT', set.ident = TRUE)

# Visualize the distribution of cell cycle markers across
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_CellCycle_RidgePlot.pdf"), width = 20, height = 20)
RidgePlot(seurat.obj.big, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
dev.off()

# calculate the difference between G2M and S pahse scores and regress it out
seurat.obj.big$CC.Difference <- seurat.obj.big$S.Score - seurat.obj.big$G2M.Score

# view cell cycle scores and phase assignments
head(seurat.obj.big[[]])

# normalise again but this time including also the cell cycle scores
seurat.obj.big <- SCTransform(seurat.obj.big, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'CC.Difference'))

# When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
seurat.obj.big <- RunPCA(seurat.obj.big, features = c(s.genes, g2m.genes))

# save with regressed out cell cycle difference
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_CellCycle_Dimplot.pdf"), width = 20, height = 20)
print(DimPlot(seurat.obj.big, reduction = "pca"))
dev.off()

# save intermediate object
saveRDS(seurat.obj.big, file = paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_SCtransform_seurat_obj.rds"))

#############################################################################################################################

## DIMENSIONALITY REDUCTION
#' With both assays in place we can now run a pca to reduce the dimensionality of the data to our essential components. Each of the analysis
#' steps will be implemented on both assays, so that we can decide which one to use at the very end. Within the next part of this document,
#' we will evaluate how many prinicipal components of the dataset accurately explain most of the variance we observe in our data. To do so, we 
#' will use `Seurat's` visualization functions and manually investigate how many PCs are meaningful for downstream analysis. One of the functions,
#' `DimHeatmap` allows for easy exploration of the primary sources of heterogeneity in a dataset cells and features are ordered according to there PCA scores.
#' Setting the argument 'cells' shows the 'extreme' cells on both ends of the spectrum. In the end, we will also base our decision on the examination
#' of an ElbowPlot, where PCs are ranked based on the percentage of variance they explain. The number of PCs will then be used for the classification
#' of clusters and nearest neighbors downstream.

## SCT ASSAY
# run principal component analysis
seurat.obj.big <- RunPCA(seurat.obj.big, npcs = 100, features = seurat.obj.big@assays$SCT@var.features)

# perform Independent component analysis as well
seurat.obj.big <- RunICA(seurat.obj.big, features = seurat.obj.big@assays$SCT@var.features)

# examine and visualize the PCA results in different ways
print(seurat.obj.big[["pca"]], dims = 1:5, nfeatures = 5)

# show PCA components with variable features
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_PCA_VizualizationPlots.pdf"), width = 20, height = 20)
print(VizDimLoadings(seurat.obj.big, dims = 1:2, reduction = "pca"))
print(DimPlot(seurat.obj.big, reduction = "pca"))
dev.off()

# make dim heatmaps
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_PCA_heatmaps.pdf"), width = 20, height = 20)
print(DimHeatmap(seurat.obj.big, dims = 1, cells = 500, balanced = T))
print(DimHeatmap(seurat.obj.big, dims = 1:20, cells = 500, balanced = T))
dev.off()

# show ICA components with variable features
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_ICA_VizualizationPlots.pdf"), width = 20, height = 20)
print(VizDimLoadings(seurat.obj.big, dims = 1:3, reduction = "ica"))
dev.off()

# Heuristic method - Elbow plot - ranking of pcs on the percentage of variance explained
# looking for an elbow, which should describe the number of pcs that describe the majority of true signal
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_ElbowPlot.pdf"), width = 20, height = 20)
print(ElbowPlot(seurat.obj.big, ndims = 100))
dev.off()

#############################################################################################################################

## CELL CLUSTERING AND EMBEDDINGS
#' With the dimensioality reduction in place, we can apply graph-based clustering method where we construct a KNN graph based
#' on the euclidean distance in PCA space. We further use modularity optimization techniques such as the louvain algorithm. In brief,
#' the higher the resolution we set in this clustering, the greater the number of clusters generally is. A resolution of 0.4 - 1.2
#' typically give good results for 3k cells. Generally, the goal of cell clustering and embedding is to learn the underlying manifold
#' of the data in order to place similar cells together in a low-dimensional space. 

## SCT
# find NNs
seurat.obj.big<- FindNeighbors(seurat.obj.big, dims = 1:75, assay = "SCT", reduction = "pca")

# find clusters
seurat.obj.big <- FindClusters(seurat.obj.big, resolution = 0.7)

# run umap
seurat.obj.big <- RunUMAP(seurat.obj.big, dims = 1:75, assay = "SCT", reduction = "pca")

# plot individual clusters
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_UMAPplot.pdf"))
print(DimPlot(seurat.obj.big, reduction = "umap", label = T))
dev.off()

# run tsne
seurat.obj.big <- RunTSNE(seurat.obj.big, reduction = "pca", dims = 1:75)

# plot tsne
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_TSNEplot.pdf"))
print(DimPlot(seurat.obj.big, reduction = "tsne", label = T))
dev.off()

# save intermediate object
saveRDS(seurat.obj.big, file = paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_DimRed_seurat_obj.rds"))

#############################################################################################################################

## DOUBLET DETECTION WITH DOUBLETFINDER
#' Uses a fully-processed seurat object (after tsne has been run) needs PCs (range 1:10 for example).
#' pN - number of artifical generated doublets (default 25%)
#' pK - definition PC neighborhodd size
#' nExp - defines threshold for doublet/singlet predictions
#' pK Identification with sctransform used (no ground-truth - can be run with ground-truth as well)
#' Further information: https://github.com/chris-mcginnis-ucsf/DoubletFinder

# follow the DoubletFinder vignette for detecting and removing doublets
sweep.res.list.seurat.obj.big <- paramSweep_v3(seurat.obj.big, PCs = 1:75, sct = T)
sweep.stats_seurat.obj.big <- summarizeSweep(sweep.res.list.seurat.obj.big, GT = FALSE)

# plot Mean-variance normalized bimodality coefficient (bcmvn) 
# ground-truth-agnostic metric that coincides with pK
bcmvn.seurat.obj.big <- find.pK(sweep.stats_seurat.obj.big)

# plot pk value
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_pKvalue_plot.pdf"))
pK=as.numeric(as.character(bcmvn.seurat.obj.big$pK))
BCmetric=bcmvn.seurat.obj.big$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
dev.off()

# get annotatios for clustering
annotations <- seurat.obj.big@meta.data$seurat_clusters

## Homotypic Doublet Proportion Estimate
homotypic.prop <- modelHomotypic(annotations) 
nExp_poi <- round(0.075*length(seurat.obj.big@meta.data$barcode))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies
seurat.obj.big <- doubletFinder_v3(seurat.obj.big, PCs = 1:75, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
seurat.obj.big <- doubletFinder_v3(seurat.obj.big, PCs = 1:75, pN = 0.25, pK = pK_choose, nExp = nExp_poi.adj, reuse.pANN = paste0("pANN_0.25_",pK_choose, "_", nExp_poi), sct = TRUE)

pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_Doublet_plot.pdf"))
print(DimPlot(seurat.obj.big, reduction = "umap", group.by = paste0("DF.classifications_0.25_",pK_choose, "_", nExp_poi.adj)))
dev.off()

# assign a new column name to the classification column
colnames(seurat.obj.big@meta.data)[ncol(seurat.obj.big@meta.data)] <- "Doublet_Score"

# remove doublets from seurat object
seurat.obj.big <- subset(seurat.obj.big, subset = Doublet_Score != "Doublet")

# save interemediate objcet
saveRDS(seurat.obj.big, file = paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_seurat_obj.rds"))

#############################################################################################################################

## SUBSETTING THE DATASET TO ONLY EPITHELIAL CELLS
#' The Zhang dataset comprises not only epithelial cells, but also stromal and immune cells. However, our organoids do not 
#' comprise these other cells, which is why we need to exclude them for further comparison of the features between 
#' the organoid and tissue epithelial cells. Thus, we will subset the analysed Zhang dataset to epithelial cells only, then
#' analysing and inferring cell types from it.

# change the cluster identities from RTAM to SCT
Idents(seurat.obj.big) <- "SCT_snn_res.0.7"

## EPITHELIAL AND NON-EPITHELIAL
# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_EpithelialMarkers_Umap_Featureplot.pdf"), width = 20, height = 20)
FeaturePlot(seurat.obj.big, features = c("EPCAM", "KRT18", "MUC1", "KRT19", "CDH1", "CLDN4"))
FeaturePlot(seurat.obj.big, features = c("CD4", "VIM", "ACTA2", "PTPRC"))
dev.off()

pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_EpithelialMarkers_expressionPlot.pdf"), width = 20, height = 20)
print(VlnPlot(object = seurat.obj.big, features = c("EPCAM", "KRT18", "MUC1", "KRT19", "CDH1", "CLDN4")))
print(VlnPlot(object = seurat.obj.big, features = c("CD4", "VIM", "ACTA2", "PTPRC", "VWF")))
dev.off()

# set idents to epithelial or non-epithelial
new.cluster.ids <- c("Non-Epithelial", "Epithelial", "Non-Epithelial", "Non-Epithelial", "Non-Epithelial", "Non-Epithelial", "Non-Epithelial", "Epithelial", "Epithelial", 
                     "Non-Epithelial", "Non-Epithelial", "Non-Epithelial", "Non-Epithelial", "Epithelial", "Non-Epithelial", "Non-Epithelial", "Epithelial", "Epithelial",
                     "Non-Epithelial", "Non-Epithelial", "Non-Epithelial", "Non-Epithelial", "Non-Epithelial", "Non-Epithelial", "Non-Epithelial", "Epithelial", "Epithelial", 
                     "Epithelial", "Epithelial", "Epithelial")

names(new.cluster.ids) <- levels(seurat.obj.big)
seurat.obj.big <- RenameIdents(seurat.obj.big, new.cluster.ids)

# subset complete seurat object to epithelial cells
seurat.obj.big$EpithelialCelltype <- Idents(seurat.obj.big)
seurat.obj.big <- subset(seurat.obj.big, subset = EpithelialCelltype == "Epithelial")

# perform PCA on variable features
seurat.obj.big <- RunPCA(seurat.obj.big, npcs = 100, features = seurat.obj.big@assays$SCT@var.features)

# alternative heuristic method - Elbow plot - ranking of pcs on the percentage of variance explained
# looking for an elbow, which should describe the number of pcs that describe the majority of true signal
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_Epithelial_ElbowPlot.pdf"), width = 20, height = 20)
print(ElbowPlot(seurat.obj.big, ndims = 100))
dev.off()

# first, construct a KNN graph based on the euclidean distance in PCA space
seurat.obj.big <- FindNeighbors(seurat.obj.big, dims = 1:50)

# apply modulairty optimization techniques 
seurat.obj.big <- FindClusters(seurat.obj.big, resolution = 0.5)

# run umap
seurat.obj.big <- RunUMAP(seurat.obj.big, dims = 1:50, assay = "SCT", reduction = "pca")

# plot individual clusters
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_Epithelial_UMAPplot.pdf"))
print(DimPlot(seurat.obj.big, reduction = "umap", label = T))
dev.off()

# run tsne
seurat.obj.big <- RunTSNE(seurat.obj.big, reduction = "pca", dims = 1:50)

# plot tsne
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_Epithelial_TSNEplot.pdf"))
print(DimPlot(seurat.obj.big, reduction = "tsne", label = T))
dev.off()

# save epithelial seurat object
saveRDS(seurat.obj.big, file = paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_Epithelial_seurat_obj.rds"))

#############################################################################################################################

## DIFFERENTIAL GENE EXPRESSION ANALYSIS
#' In this last part, we will create some plots and investigate distinct clusters for the expression of different genes, most
#' of which have been reported in the literature. We will use this manual investigation, to manually, but only roughly, 
#' assign possible identities to clusters. 

# # read seurat objects
# seurat.obj.big <- readRDS(file = paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_Epithelial_seurat_obj.rds"))

## CELL MARKERS FROM ZHANG ET AL.
# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_featurePlot.pdf"), width = 20, height = 20)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("MUC5AC", "TFF1", "CD79A", "CD19", "OLFM4", "LGR5", "SOX2","CCKBR", "FABP1", "FABP2", "CA1", "VIL1", "MUC6", "TFF2", "CHGA", "TAC1", "TPH1", "CHGB"))
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("CDK1", "MKI67", "CEACAM5", "CEACAM6", "PGC", "CXCL3", "IL8", "COL1A2", "LUM", "DCN", "PDPN", "FAP", "COL3A1", "COL6A1", "VWF", "ENG", "MCAM"))
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("SPINK4", "TFF3", "MUC2", "ITLN1", "CD14", "CD68", "CSF1R", "MYL2", "ACTA2", "CD14", "CD68", "CSF1R", "MYL2"))
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("ACTA2", "TPSAB1", "TPSB2", "PGA3", "PGA4", "LIPF", "CD2", "CD3D", "CD3E", "CD3G", "ATP4A", "ATP4B", "GAST", "GHRL", "SST"))
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("OLFM4", "PHLDA1", "LEFTY1")) # stem cells
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("CEACAM6", "BAX", "CCND2")) # cancer cells
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("CEACAM5", "FABP1", "CDH17")) # non-specific cancer cells - also in enterocytes
dev.off()

## EPITHELIAL AND NON-EPITHELIAL
# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_EpithelialMarkers_Umap_Featureplot.pdf"), width = 20, height = 20)
FeaturePlot(seurat.obj.big, features = c("EPCAM", "KRT18", "MUC1", "KRT19", "CDH1", "CLDN4"))
FeaturePlot(seurat.obj.big, features = c("CD4", "VIM", "ACTA2", "PTPRC"))
dev.off()

pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_EpithelialMarkers_expressionPlot.pdf"), width = 20, height = 20)
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("EPCAM", "KRT18", "MUC1", "KRT19", "CDH1", "CLDN4")))
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("CD4", "VIM", "ACTA2", "PTPRC")))
dev.off()

# plot probability distributions across clusters for the top5 genes indicated in the Supplementary data from Zhang et al.
# However, markers shared by different cell types were removed from both.
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_exp_violin_plot.pdf"), width = 20, height = 20)
# check for PMCs
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("GKN1", "GKN2", "MUC5AC", "TFF1", "DPCR1")) + labs(title = "PMCs"))
# check for MSCs
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("OLFM4", "RPL7", "CLDN4","TSPAN8", "REG1A")) + labs(title = "MSCs"))
# check for enterocytes
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("FABP1", "FABP2", "RBP2", "ANPEP", "APOA4")) + labs(title = "Enterocytes"))
# check for GMC
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("PRR4", "C6orf58","MUC6", "TFF2", "LTF")) + labs(title = "GMCs"))
# check for enteroendocrine
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("CHGA", "PCSK1N", "SCG5", "CHGB", "TPH1")) + labs(title = "Enteroendocrine")) 
# check for PCs
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("TOP2A", "MKI67", "UBE2C", "HMGB2", "PTTG1")) + labs(title = "PCs"))
# check for cancer cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("REG4","CLDN7","KRT18", "LGALS3", "CEACAM6")) + labs(title = "Cancer Cells"))
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("CLDN3", "CST1", "MUC3A", "CLDN4", "PI3", "UBD")) + labs(title = "Cancer Cells"))
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("CDH17", "PRAP1", "UBE2C", "CCL20", "LCN2", "SERPINB5")) + labs(title = "Cancer Cells"))
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("RRM2", "MYBL2", "MMP7", "TPX2", "MISP", "TMPRSS4")) + labs(title = "Cancer Cells"))
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("RMRP", "CLDN1", "GPRC5A", "CLRN3", "CXCL1", "MSLN")) + labs(title = "Cancer Cells"))
# check for neck like cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("MIA", "CXCL3", "CXCL2", "CXCL17", "CLU")) + labs(title = "Neck-like cells"))
# check for Goblet cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("SPINK4", "TFF3", "MUC2", "ITLN1", "ZG16")) + labs(title = "Goblet Cells"))
# check for Chief cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("PGA3", "PGA4", "LIPF", "CHIA", "PDIA2")) + labs(title = "Chief Cells")) 
# check for parietal cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("ATP4A", "ATP4B", "VEGFB")) + labs(title = "Parietal cells"))
# check for endocrine cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("GAST", "GHRL", "SST")) + labs(title = "Endocrine cells"))
##GAO DATASET
# check for Chief Stem Cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("GIF", "LRIG1", "PROCR")) + labs(title = "Chief Stem cells"))
# check for Enteroendocrine
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("SST", "PYY", "PROX1", "TPH1", "REG4", "NEUROD1", "GIP")) + labs(title = "Enteroendocrine cells"))
# check for Enterocytes
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("CA2", "KRT20", "ALPI", "TREH", "LCT", "MME", "CDH1", "VIL1", "AQP8", "SI", "CDX2")) + labs(title = "Enterocytes"))
# check for paneth cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("NOTUM", "NUPR1", "PLA2G2A")) + labs(title = "Paneth Cells"))
# check for tuft cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("DCLKK1", "TRPM5", "PTGS1", "RGS13")) + labs(title = "Tuft Cells"))
# check for Stem cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("PROM1", "MEX3A", "TERT", "LGR5", "SOX9", "CD24", "ALCAM", "PROCR")) + labs(title = "Stem Cells"))
# check for PROCRhigh progenitor cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("PPIH", "SGK1", "GPT2")) + labs(title = "PROCR High Progenitor Cells"))
# check for Neck progenitor cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("TFF2", "ECE1", "PROM1")) + labs(title = "Neck-Progenitor Cells"))
# check for HES1high progenitor cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("LINC01578", "PLXDC1")) + labs(title = "HES1 High Progenitor Cells"))
# check for parietal progenitor cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("ATP4A", "ATP4B", "CFH", "FGD5")) + labs(title = "Parietal Progenitor Cells"))
# check for Pit progenitor cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("LRIG1", "AXIN2", "CD44", "ACTC1")) + labs(title = "Pit Progenitor Cells"))
# check for Pit Cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("ID1", "MBD1", "GREM", "RSPO2")) + labs(title = "Pit Cells"))
dev.off()

###########################################################################
#               ANALYSE THE GENE EXPRESSION PER CLUSTER
###########################################################################

# find markers for every cluster compared to all remaining cells, report only the positive ones
seurat.obj.big.markers <- FindAllMarkers(seurat.obj.big, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5_seurat.obj.big.markers <- seurat.obj.big.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
top10_seurat.obj.big.markers <- seurat.obj.big.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

# write table with the DGEs per cluster
write.table(seurat.obj.big.markers, paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_cluster_markers.txt"), sep = "\t", quote = F, col.names = T, row.names = F)

## TOP MARKERS
# top 5 markers (or all markers if less than 10) for each cluster.
pdf(paste(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_ExpressionHeatmap", ".pdf", sep = ""), width = 10.5, height = 9)
print(DoHeatmap(seurat.obj.big, features = top5_seurat.obj.big.markers$gene, size = 3, angle = 45, disp.min = -3, disp.max = 3) + scale_fill_gradientn(colors = c("blue", "lightgrey", "red")))
dev.off()

pdf(paste(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_dataSlot_ExpressionHeatmap", ".pdf", sep = ""), width = 10.5, height = 9)
print(DoHeatmap(seurat.obj.big, assay = "SCT", slot = "data", features = top5_seurat.obj.big.markers$gene, size = 3, angle = 45, disp.max = 3) + scale_fill_gradientn(colors = c("white", "blue")))
dev.off()

############################################################################
##    ANNOTATE THE EL ORGANOID DATA USING THE LITERATURE MARKER LIST
############################################################################

# Idents(seurat.obj.big) <- "seurat_clusters"

# use canonical markers to match the unbiased clustering to known cell types for Zhang dataset
new.cluster.ids <- c("Malignant Cells", "Enterocytes", "Chief Cells", "PMCs", "Goblet Cells", "Malignant Cells", "Unknown", "Unknown",
                     "Enterocytes", "MSCs", "Chief Cells", "Enteroendocrine", "GMCs", "Unknown")

# change cluster ids for Cell names
names(new.cluster.ids) <- levels(seurat.obj.big)
seurat.obj.big <- RenameIdents(seurat.obj.big, new.cluster.ids)

# subset complete seurat object to epithelial cells
seurat.obj.big$Celltype <- Idents(seurat.obj.big)
seurat.obj.big <- subset(seurat.obj.big, subset = Celltype != "Unknown")

# perform PCA on variable features
seurat.obj.big <- RunPCA(seurat.obj.big, npcs = 100, features = seurat.obj.big@assays$SCT@var.features)

# alternative heuristic method - Elbow plot - ranking of pcs on the percentage of variance explained
# looking for an elbow, which should describe the number of pcs that describe the majority of true signal
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_Epithelial_Celltype_ElbowPlot.pdf"), width = 20, height = 20)
print(ElbowPlot(seurat.obj.big, ndims = 100))
dev.off()

# first, construct a KNN graph based on the euclidean distance in PCA space
seurat.obj.big <- FindNeighbors(seurat.obj.big, dims = 1:50)

# apply modulairty optimization techniques 
seurat.obj.big <- FindClusters(seurat.obj.big, resolution = 0.5)

# run umap
seurat.obj.big <- RunUMAP(seurat.obj.big, dims = 1:50, assay = "SCT", reduction = "pca")

# plot individual clusters
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_Epithelial_Celltype_UMAPplot.pdf"))
print(DimPlot(seurat.obj.big, reduction = "umap", label = T))
dev.off()

# run tsne
seurat.obj.big <- RunTSNE(seurat.obj.big, reduction = "pca", dims = 1:50)

# plot tsne
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_Epithelial_Celltype_TSNEplot.pdf"))
print(DimPlot(seurat.obj.big, reduction = "tsne", label = T))
dev.off()

# save epithelial seurat object
saveRDS(seurat.obj.big, file = paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_Epithelial_seurat_obj.rds"))

#############################################################################################################################

## DIFFERENTIAL GENE EXPRESSION ANALYSIS
#' In this last part, we will create some plots and investigate distinct clusters for the expression of different genes, most
#' of which have been reported in the literature. We will use this manual investigation, to manually, but only roughly, 
#' assign possible identities to clusters. 

## CELL MARKERS FROM ZHANG ET AL.
# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_featurePlot.pdf"), width = 20, height = 20)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("MUC5AC", "TFF1", "CD79A", "CD19", "OLFM4", "LGR5", "SOX2","CCKBR", "FABP1", "FABP2", "CA1", "VIL1", "MUC6", "TFF2", "CHGA", "TAC1", "TPH1", "CHGB"))
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("CDK1", "MKI67", "CEACAM5", "CEACAM6", "PGC", "CXCL3", "IL8", "COL1A2", "LUM", "DCN", "PDPN", "FAP", "COL3A1", "COL6A1", "VWF", "ENG", "MCAM"))
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("SPINK4", "TFF3", "MUC2", "ITLN1", "CD14", "CD68", "CSF1R", "MYL2", "ACTA2", "CD14", "CD68", "CSF1R", "MYL2"))
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("ACTA2", "TPSAB1", "TPSB2", "PGA3", "PGA4", "LIPF", "CD2", "CD3D", "CD3E", "CD3G", "ATP4A", "ATP4B", "GAST", "GHRL", "SST"))
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("OLFM4", "PHLDA1", "LEFTY1")) # stem cells
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("CEACAM6", "BAX", "CCND2")) # cancer cells
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("CEACAM5", "FABP1", "CDH17")) # non-specific cancer cells - also in enterocytes
dev.off()

## EPITHELIAL AND NON-EPITHELIAL
# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_EpithelialMarkers_Umap_Featureplot.pdf"), width = 20, height = 20)
FeaturePlot(seurat.obj.big, features = c("EPCAM", "KRT18", "MUC1", "KRT19", "CDH1", "CLDN4"))
FeaturePlot(seurat.obj.big, features = c("CD4", "VIM", "ACTA2", "PTPRC"))
dev.off()

pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_EpithelialMarkers_expressionPlot.pdf"), width = 20, height = 20)
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("EPCAM", "KRT18", "MUC1", "KRT19", "CDH1", "CLDN4")))
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("CD4", "VIM", "ACTA2", "PTPRC")))
dev.off()

# plot probability distributions across clusters for the top5 genes indicated in the Supplementary data from Zhang et al.
# However, markers shared by different cell types were removed from both.
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_exp_violin_plot.pdf"), width = 20, height = 20)
# check for PMCs
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("GKN1", "GKN2", "MUC5AC", "TFF1", "DPCR1")) + labs(title = "PMCs"))
# check for MSCs
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("OLFM4", "RPL7", "CLDN4","TSPAN8", "REG1A")) + labs(title = "MSCs"))
# check for enterocytes
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("FABP1", "FABP2", "RBP2", "ANPEP", "APOA4")) + labs(title = "Enterocytes"))
# check for GMC
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("PRR4", "C6orf58","MUC6", "TFF2", "LTF")) + labs(title = "GMCs"))
# check for enteroendocrine
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("CHGA", "PCSK1N", "SCG5", "CHGB", "TPH1")) + labs(title = "Enteroendocrine")) 
# check for PCs
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("TOP2A", "MKI67", "UBE2C", "HMGB2", "PTTG1")) + labs(title = "PCs"))
# check for cancer cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("REG4","CLDN7","KRT18", "LGALS3", "CEACAM6")) + labs(title = "Cancer Cells"))
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("CLDN3", "CST1", "MUC3A", "CLDN4", "PI3", "UBD")) + labs(title = "Cancer Cells"))
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("CDH17", "PRAP1", "UBE2C", "CCL20", "LCN2", "SERPINB5")) + labs(title = "Cancer Cells"))
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("RRM2", "MYBL2", "MMP7", "TPX2", "MISP", "TMPRSS4")) + labs(title = "Cancer Cells"))
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("RMRP", "CLDN1", "GPRC5A", "CLRN3", "CXCL1", "MSLN")) + labs(title = "Cancer Cells"))
# check for neck like cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("MIA", "CXCL3", "CXCL2", "CXCL17", "CLU")) + labs(title = "Neck-like cells"))
# check for Goblet cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("SPINK4", "TFF3", "MUC2", "ITLN1", "ZG16")) + labs(title = "Goblet Cells"))
# check for Chief cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("PGA3", "PGA4", "LIPF", "CHIA", "PDIA2")) + labs(title = "Chief Cells")) 
# check for parietal cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("ATP4A", "ATP4B", "VEGFB")) + labs(title = "Parietal cells"))
# check for endocrine cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("GAST", "GHRL", "SST")) + labs(title = "Endocrine cells"))
##GAO DATASET
# check for Chief Stem Cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("GIF", "LRIG1", "PROCR")) + labs(title = "Chief Stem cells"))
# check for Enteroendocrine
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("SST", "PYY", "PROX1", "TPH1", "REG4", "NEUROD1", "GIP")) + labs(title = "Enteroendocrine cells"))
# check for Enterocytes
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("CA2", "KRT20", "ALPI", "TREH", "LCT", "MME", "CDH1", "VIL1", "AQP8", "SI", "CDX2")) + labs(title = "Enterocytes"))
# check for paneth cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("NOTUM", "NUPR1", "PLA2G2A")) + labs(title = "Paneth Cells"))
# check for tuft cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("DCLKK1", "TRPM5", "PTGS1", "RGS13")) + labs(title = "Tuft Cells"))
# check for Stem cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("PROM1", "MEX3A", "TERT", "LGR5", "SOX9", "CD24", "ALCAM", "PROCR")) + labs(title = "Stem Cells"))
# check for PROCRhigh progenitor cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("PPIH", "SGK1", "GPT2")) + labs(title = "PROCR High Progenitor Cells"))
# check for Neck progenitor cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("TFF2", "ECE1", "PROM1")) + labs(title = "Neck-Progenitor Cells"))
# check for HES1high progenitor cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("LINC01578", "PLXDC1")) + labs(title = "HES1 High Progenitor Cells"))
# check for parietal progenitor cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("ATP4A", "ATP4B", "CFH", "FGD5")) + labs(title = "Parietal Progenitor Cells"))
# check for Pit progenitor cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("LRIG1", "AXIN2", "CD44", "ACTC1")) + labs(title = "Pit Progenitor Cells"))
# check for Pit Cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("ID1", "MBD1", "GREM", "RSPO2")) + labs(title = "Pit Cells"))
dev.off()

###########################################################################
#               ANALYSE THE GENE EXPRESSION PER CLUSTER
###########################################################################

# find markers for every cluster compared to all remaining cells, report only the positive ones
seurat.obj.big.markers <- FindAllMarkers(seurat.obj.big, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5_seurat.obj.big.markers <- seurat.obj.big.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
top10_seurat.obj.big.markers <- seurat.obj.big.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

# write table with the DGEs per cluster
write.table(seurat.obj.big.markers, paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_cluster_markers.txt"), sep = "\t", quote = F, col.names = T, row.names = F)

## TOP MARKERS
# top 5 markers (or all markers if less than 10) for each cluster.
pdf(paste(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_ExpressionHeatmap", ".pdf", sep = ""), width = 10.5, height = 9)
print(DoHeatmap(seurat.obj.big, features = top5_seurat.obj.big.markers$gene, size = 3, angle = 45, disp.min = -3, disp.max = 3) + scale_fill_gradientn(colors = c("blue", "lightgrey", "red")))
dev.off()

pdf(paste(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_dataSlot_ExpressionHeatmap", ".pdf", sep = ""), width = 10.5, height = 9)
print(DoHeatmap(seurat.obj.big, assay = "SCT", slot = "data", features = top5_seurat.obj.big.markers$gene, size = 3, angle = 45, disp.max = 3) + scale_fill_gradientn(colors = c("white", "blue")))
dev.off()

############################################################################
##    ANNOTATE THE EL ORGANOID DATA USING THE LITERATURE MARKER LIST
############################################################################

# Idents(seurat.obj.big) <- "seurat_clusters"

# use canonical markers to match the unbiased clustering to known cell types for Zhang dataset
new.cluster.ids <- c("Malignant Cells", "Enterocytes", "PMCs", "Chief Cells", "Malignant Cells", "Goblet Cells",
                     "Enterocytes", "MSCs", "Parietal Cells", "Chief Cells","Enteroendocrine", "GMCs", "Goblet Cells")

# change cluster ids for Cell names
names(new.cluster.ids) <- levels(seurat.obj.big)
seurat.obj.big <- RenameIdents(seurat.obj.big, new.cluster.ids)
seurat.obj.big$Celltype <- Idents(seurat.obj.big)

###########################################################################
#             ANALYSE THE GENE EXPRESSION PER CLONE
###########################################################################

# get markers per celltype
celltype.seurat.obj.big.markers <- FindAllMarkers(seurat.obj.big, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10_celltype_seurat.obj.big.markers <- celltype.seurat.obj.big.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

# write table with the DGEs per cluster
write.table(celltype.seurat.obj.big.markers, paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_Celltype_markers.txt"), sep = "\t", quote = F, col.names = T, row.names = F)

## TOP MARKERS
# top 3 markers (or all markers if less than 10) for each cluster.
pdf(paste(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_Celltype_ExpressionHeatmap", ".pdf", sep = ""), width = 10.5, height = 9)
print(DoHeatmap(seurat.obj.big, features = top10_celltype_seurat.obj.big.markers$gene, size = 3, angle = 45, disp.min = -3, disp.max = 3) + scale_fill_gradientn(colors = c("blue", "lightgrey", "red")))
dev.off()

# top 3 markers (or all markers if less than 10) for each cluster.
pdf(paste(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_dataSlot_Celltype_ExpressionHeatmap", ".pdf", sep = ""), width = 10.5, height = 9)
print(DoHeatmap(seurat.obj.big, assay = "SCT", slot = "data", features = top10_celltype_seurat.obj.big.markers$gene, size = 3, angle = 45, disp.max = 3) + scale_fill_gradientn(colors = c("white", "blue")))
dev.off()

############################################################################
##        MAKE PLOTS FOR THE UMAP EMBEDDING FOR EL ORGANOIDS
############################################################################

## CLONE ID
plot1 <- DimPlot(seurat.obj.big, reduction = "umap", group.by = "seurat_clusters", pt.size = 1.5, label = T, label.size = 6) +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=20, face = "bold")) + 
  theme(legend.text = element_text(colour="black", size=20, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  guides(color=guide_legend("Clusters", override.aes = list(size = 6)))

# save single plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_Clusters.pdf"), width = 6, height = 5)
print(plot1 + NoLegend()) 
dev.off()

plot2 <- DimPlot(seurat.obj.big, reduction = "umap", group.by = "timepoint", pt.size = 1.5) +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=40, face = "bold")) + 
  theme(legend.text = element_text(colour="black", size=40, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  scale_color_manual(labels = c("Early", "Late", "WT"), values = c("#EA9F37", "#782867", "#4F9E4C")) +
  guides(color=guide_legend("Timepoint", override.aes = list(size = 8)))

# save single plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_Timepoint.pdf"), width = 6, height = 5)
print(plot2 + NoLegend()) 
dev.off()

plot3 <- DimPlot(seurat.obj.big, reduction = "umap", pt.size = 1.5, label = T, label.size = 3) +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=20, face = "bold")) + 
  theme(legend.text = element_text(colour="black", size=20, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  guides(color=guide_legend("Cell Types", override.aes = list(size = 6)))

# save single plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_CellTypes.pdf"), width = 6, height = 5)
print(plot3 + NoLegend()) 
dev.off()

plotlist <- list(plot1, plot2, plot3)

# plot the list
nCol <- 3
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_Comparison_plot.pdf"), width = 60, height = 20)
cowplot::plot_grid(plotlist = plotlist, ncol = nCol)
dev.off()

# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_GenesOfInterest1_FeaturePlot.pdf"), width = 11, height = 9)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("MUC5AC", "TFF1", "CEACAM6", "PGC"), pt.size = 1, ncol = 2)
dev.off()

# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_GenesOfInterest2_FeaturePlot.pdf"), width = 11, height = 9)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("WFDC2", "MUC5B", "REG4", "OLFM4"), pt.size = 1, ncol = 2)
dev.off()

# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_GenesOfInterest3_FeaturePlot.pdf"), width = 11, height = 9)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("TSPAN8", "VIL1", "AKR1C1", "LYZ"), pt.size = 1, ncol = 2)
dev.off()

# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_GenesOfInterest4_FeaturePlot.pdf"), width = 11, height = 9)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("DMBT1", "FABP1", "GSTA1", "PSCA"), pt.size = 1, ncol = 2)
dev.off()

# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_GeneHeterogeneity_FeaturePlot.pdf"), width = 11, height = 9)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("MUC5AC", "TFF1", "CEACAM6", "WFDC2", "REG4", "DMBT1", "FABP1", "GSTA1", "LYZ"), pt.size = 1, ncol = 3)
dev.off()


# save final cell type annotated seurat_obj
saveRDS(seurat.obj.big, file = paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_Epithelial_CellType_seurat_obj.rds"))
