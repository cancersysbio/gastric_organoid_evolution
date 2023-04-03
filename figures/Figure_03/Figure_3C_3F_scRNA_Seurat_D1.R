#############################################################################################################################
##                                                                                                                      
##  ANALYSIS DONOR 1 EARLY MID LATE SAMPLES WITH SEURAT
##                                                                                                                      
##  Date: 14 MAY 2021                                                                                                                    
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
list.of.packages <- c("tidyverse", "patchwork", "Seurat", "Matrix", "biomaRt", "scater", "DoubletFinder", "viridis",
                      "ggsignif", "SeuratWrappers", "harmony", "future")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages, repos = "http://cran.us.r-project.org")

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

# Parallelization (t=16)
# future::plan("multiprocess", workers = 16)

# set up functions which are used 
`%notin%` <- Negate(`%in%`)

#############################################################################################################################

## READ IN YOUR DATA OF INTEREST
#' Here, we read in all the 10x Genomics matrices, that belong to the set of samples we want to merge. 
#' Merging data is only recommended if you expect lower batch effects than in for instance the analysis of datasets from distinct
#' laboratories and technologies. Here, we assume that we have minimal batch effects in our data from the same set of experiments
#' and technologies. For the dataset from Zhang et al., 2019, Cell we use the raw count matrices from pre-processed Seurat objects. 

## EARLY LATE ORGANOIDS
c.dir <- "/labs/ccurtis2/mjprzy/infercnv_gastric"

# define output directory
o.dir <- "/labs/ccurtis2/mjprzy/scRNA_analysis/hashECB_data_freeze"
setwd(o.dir)

# define a sample id for the merged object
sample.tmp <- "EML_Seq_D1"

# create sample-specific output directory with a "plots" subfolder
dir.create(paste0(o.dir, "/" , sample.tmp))
dir.create(paste0(o.dir, "/" , sample.tmp, "/plots"))

# list to store seurat.obj.bigects in
seurat.obj.big.list <- list()

# set sample ids for the seurat objects
sample.id.list <- c("Seq27", "Seq29", "Seq30")

i <- 1
# iterate over each matrix, subset to cells of interest and create seurat objects
for (i in 1:length(sample.id.list)){
  
  print(sample.id.list[i])
  load(file = paste0(sample.tmp, "/sc.10x.counts_", sample.id.list[i], ".RData"))

  # read
  metadata <- read.table(paste0("EML_Seq_D1/", sample.id.list[i], "_metadata.txt"), header = T)
  
  # subset matrices to cells which are in metadata
  c.mtx <- c.mtx[,colnames(c.mtx) %in% metadata$Cell_Barcode]
  c.mtx <- c.mtx[,unique(colnames(c.mtx))]
  
  # check dims
  print(dim(c.mtx))
  
  # create new seurat object
  seurat.obj.big <- CreateSeuratObject(counts = c.mtx, 
                                       project = sample.id.list[i], 
                                       min.cells = 3, 
                                       min.features = 200)
  
  # add Cell_barcode column
  seurat.obj.big$Cell_Barcode <- rownames(seurat.obj.big@meta.data)
  
  # merge with the Hash metadata
  seurat.obj.big@meta.data <- merge(seurat.obj.big@meta.data, metadata, by = "Cell_Barcode")
  
  # store in list
  seurat.obj.big.list[[i]] <- seurat.obj.big
  
}

# combined the seurat objects together
seurat.obj.big <- merge(seurat.obj.big.list[[1]], y = seurat.obj.big.list[c(2:length(seurat.obj.big.list))], 
                        add.cell.ids = sample.id.list, project = "EML_Seq", merge.data = T)

# save identity
seurat.obj.big$sample_ident <- seurat.obj.big$orig.ident

# replace the orig.ident column with HashTag
seurat.obj.big$orig.ident <- seurat.obj.big$HashTag

# check colnames
head(colnames(seurat.obj.big))
tail(colnames(seurat.obj.big))

# check amount of cells
table(seurat.obj.big$orig.ident)

# check the cell contributions before QC
pre.data <- as.data.frame(table(seurat.obj.big$orig.ident))
colnames(pre.data) <- c("Sample", "Cells_pre")

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

# check the cell contributions before QC
post.data <- as.data.frame(table(seurat.obj.big$orig.ident))
colnames(post.data) <- c("Sample", "Cells_post")

# pre.data
adapted.pre.data <- pre.data
adapted.pre.data$Cells_pre <- adapted.pre.data$Cells_pre - post.data$Cells_post

# visualize the number of cells pre and post QC
complete.cell.data <- merge(adapted.pre.data, post.data, by = "Sample")
complete.cell.data.melt <- reshape2::melt(complete.cell.data)

# Small multiple
ggplot(complete.cell.data.melt, aes(y=value, x=Sample, fill = variable)) + 
  geom_bar(position="stack", stat="identity") +
  xlab("") + theme_classic() + ylab("# mtDNA mutations") +
  scale_fill_brewer(palette = "Greens") +
  theme(axis.text = element_text(colour = "black", size = 10, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 45, hjust = 1),
        axis.title = element_text(colour = "black", size = 12, face = "bold" ),
        plot.title = element_text(colour = "black", size = 12, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 10, face = "bold",),
        legend.text = element_text(colour="black", size=10, face="bold"), 
        legend.position = "top")
ggsave(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_pre_post_cells.pdf"), width = 10, height = 5)

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
seurat.obj.big <- SCTransform(seurat.obj.big, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c('sample_ident', 'percent.mt', 'nFeature_RNA', 'nCount_RNA'))

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
seurat.obj.big <- SCTransform(seurat.obj.big, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c('sample_ident', 'percent.mt', 'nFeature_RNA', 'nCount_RNA', 'CC.Difference'))

# When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
seurat.obj.big <- RunPCA(seurat.obj.big, features = c(s.genes, g2m.genes))

# save with regressed out cell cycle difference
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_CellCycle_Dimplot.pdf"), width = 20, height = 20)
print(DimPlot(seurat.obj.big, reduction = "pca"))
dev.off()

# save intermediate object
saveRDS(seurat.obj.big, file = paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_SCtransform_Patient_seurat_obj.rds"))

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
#' 
## SCT ASSAY
# run principal component analysis
seurat.obj.big <- RunPCA(seurat.obj.big, npcs = 50, features = seurat.obj.big@assays$SCT@var.features)

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
print(ElbowPlot(seurat.obj.big, ndims = 50))
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
seurat.obj.big<- FindNeighbors(seurat.obj.big, dims = 1:15, assay = "SCT", reduction = "pca")

# find clusters
seurat.obj.big <- FindClusters(seurat.obj.big, resolution = 0.8)

# run umap
seurat.obj.big <- RunUMAP(seurat.obj.big, dims = 1:15, assay = "SCT", reduction = "pca")

# plot individual clusters
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_UMAPplot.pdf"))
print(DimPlot(seurat.obj.big, reduction = "umap", label = T))
dev.off()

# run tsne
seurat.obj.big <- RunTSNE(seurat.obj.big, reduction = "pca", dims = 1:15)

# plot tsne
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_TSNEplot.pdf"))
print(DimPlot(seurat.obj.big, reduction = "tsne", label = T))
dev.off()

# save intermediate object
saveRDS(seurat.obj.big, file = paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_DimRed_Patient_seurat_obj.rds"))

#############################################################################################################################

## DOUBLET DETECTION WITH DOUBLETFINDER
#' Uses a fully-processed seurat object (after tsne has been run) needs PCs (range 1:10 for example).
#' pN - number of artifical generated doublets (default 25%)
#' pK - definition PC neighborhodd size
#' nExp - defines threshold for doublet/singlet predictions
#' pK Identification with sctransform used (no ground-truth - can be run with ground-truth as well)
#' Further information: https://github.com/chris-mcginnis-ucsf/DoubletFinder

# follow the DoubletFinder vignette for detecting and removing doublets
sweep.res.list.seurat.obj.big <- paramSweep_v3(seurat.obj.big, PCs = 1:15, sct = T)
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
seurat.obj.big <- doubletFinder_v3(seurat.obj.big, PCs = 1:15, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
seurat.obj.big <- doubletFinder_v3(seurat.obj.big, PCs = 1:15, pN = 0.25, pK = pK_choose, nExp = nExp_poi.adj, reuse.pANN = paste0("pANN_0.25_",pK_choose, "_", nExp_poi), sct = TRUE)

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

## DIFFERENTIAL GENE EXPRESSION ANALYSIS
#' In this last part, we will create some plots and investigate distinct clusters for the expression of different genes, most
#' of which have been reported in the literature. We will use this manual investigation, to manually, but only roughly, 
#' assign possible identities to clusters. 

# set ids to clones 
seurat.obj.big$clones <- seurat.obj.big$orig.ident
unique(seurat.obj.big$clones)

## CLONE ID
plot1 <- DimPlot(seurat.obj.big, reduction = "umap", group.by = "clone_id") +
  xlab("UMAP Dimension 1") + 
  ylab("UMAP Dimension 2") +
  theme_classic() +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=25, face = "bold")) + 
  theme(legend.text = element_text(colour="black", size=25, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  guides(color=guide_legend("CLONE IDs", override.aes = list(size = 6)))

# save single plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_CloneID.pdf"), width = 6, height = 5)
print(plot1 + NoLegend() + NoAxes()) 
dev.off()

plot2 <- DimPlot(seurat.obj.big, reduction = "umap", group.by = "timepoint") +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=40, face = "bold")) + 
  theme(legend.text = element_text(colour="black", size=40, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  scale_color_manual(labels = c("Early", "Mid", "Late", "WT"), values = c("#EA9F37", "darkblue", "#782867", "#4F9E4C")) +
  guides(color=guide_legend("Timepoint", override.aes = list(size = 8)))

# save single plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_Timepoint.pdf"), width = 6, height = 5)
print(plot2 + NoLegend() + NoAxes()) 
dev.off()

## CREATE A DOTPLOT FOR THE CLONES
## WITH GASTRIC MARKERS
Idents(seurat.obj.big) <- "clones" # EL combined

levels(seurat.obj.big) <- c("WT", 
                            "C1_Early", "C1_Mid", "C1_Late",
                            "C2_Early", "C2_Mid", "C2_Late",
                            "C3_Early", "C3_Mid", "C3_Late")

plot3 <- DotPlot(seurat.obj.big, features = c("MUC5AC", "TFF1", # Pit Mucosal Genes
                                              "MUC6",  "TFF2", # Gland Mucosal Genes
                                              "MKI67", # Proliferation
                                              "PGC", "LYZ", # Neck-like Cells
                                              "OLFM4", # Mucosal Stem cells
                                              "FABP1", "VIL1", # Enterocytes
                                              "TFF3", "WFDC2", "MUC5B", "CDX2", # Goblet Cells
                                              "CEACAM5", "CEACAM6", "CLDN3", "CLDN4", "CLDN7", "REG4", "MUC13", "UBD", "AOC1", "CDH17", # GEPIA 15
                                              "TP53", "APC", "CDKN2A", "FHIT"), cols = c("lightgrey", "darkred"), scale = TRUE, col.min = 0, col.max = 3) +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=12, face = "bold")) + 
  theme(legend.text = element_text(colour="black", size=12, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=12),
        axis.text.x = element_text(face="bold", color="black", size=12, angle = 45, hjust = 1))

plot3

#############################################################################################################################
#                                     ANALYSE THE GENE EXPRESSION PER CLUSTER

# set the timepoint as identity to compare to each other
Idents(seurat.obj.big) <- "timepoint"
levels(seurat.obj.big) <- c("WT", "Early", "Mid", "Late")
unique(seurat.obj.big$timepoint)

# find markers for every cluster compared to all remaining cells, report only the positive ones
seurat.obj.big.markers <- FindAllMarkers(seurat.obj.big, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25)
seurat.obj.big.markers <- seurat.obj.big.markers[order(seurat.obj.big.markers$cluster, seurat.obj.big.markers$avg_log2FC, decreasing = c(F,T)),]
top20_seurat.obj.big.markers <- seurat.obj.big.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

# write table with the DGEs per cluster
write.table(top20_seurat.obj.big.markers, paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_top20_timepoint_markers.filtered.txt"), sep = "\t", quote = F, col.names = T, row.names = F)

## TOP MARKERS
# top 20 markers for each timepoint
pdf(paste(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_ExpressionHeatmap", ".pdf", sep = ""), width = 12, height = 9, useDingbats = T)
print(DoHeatmap(seurat.obj.big, features = top20_seurat.obj.big.markers$gene, size = 5, draw.lines = T, hjust = 0.4, angle = 0, disp.min = -3, disp.max = 3, raster = F,  group.colors = c("#4F9E4C", "#EA9F37", "darkblue", "#782867")) + 
        scale_fill_gradientn(colors = c("blue", "white", "red")) + 
        theme(text = element_text(size = 12, color = "black")))
dev.off()
