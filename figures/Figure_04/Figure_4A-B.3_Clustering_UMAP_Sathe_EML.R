#############################################################################################################################
##                                                                                                                      
##  Perform clustering and umap embedding on 10x single-cell RNA-sequencing data
##
##  Date: 16 May 2020                                                                                                                   
##  
##  Author: Moritz Przybilla
##
##  Credit: This code was designed based on code developed by Jeffrey Granja in frame of Granja*, Klemm*, Mcginnis* et al. 
##          A single cell framework for multi-omic analysis of disease identifies  malignant regulatory signatures in 
##          mixed phenotype acute leukemia (2019, https://github.com/GreenleafLab/MPAL-Single-Cell-2019). 
##
##                                                                                                                      
############################################################################################################################
# clear workspace
rm(list=ls())
set.seed(1) # set the seed of random number generator 

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("reshape2", "optparse", "BSgenome", "RColorBrewer", "ggplot2", "scales", "dendextend", "tidyverse", 
                      "Matrix", "ComplexHeatmap", "Rtsne", "robustbase", "psych", "cluster", "Matrix.utils",
                      "BSgenome.Hsapiens.UCSC.hg38", "Seurat", "GenomicRanges", "SingleCellExperiment", "matrixStats", "readr",
                      "magrittr", "edgeR", "uwot", "BiocManager", "Rcpp", "biomaRt", "dplyr", "viridis", "httr", "RColorBrewer")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

#####################################################################################
# FUNCTIONS
#####################################################################################

#' **Binarize a count matrix**
#' Use this function to convert a count matrix into a binary matrix. 
#' binarizeMat(
#'      mat = count.matrix 
#' )
#' 
#' 

#Binarize Sparse Matrix
binarizeMat <- function(mat){
  mat@x[mat@x > 0] <- 1
  mat
}

#' **Calculate Latent semantic indexing**
#' Use this function to create a GRange object that contains chr, start, end, as well as 
#' information about the GC and AT content of each window.
#' calcLSI(
#'      mat = count.matrix,
#'      nComponents = nPCs, 
#'      binarize = TRUE, 
#'      nFeatures = NULL
#' )
#' 

#LSI Adapted from fly-atac with information for re-projection analyses
calcLSI <- function(mat, nComponents = 50, binarize = TRUE, nFeatures = NULL){
  
  set.seed(1)
  
  #TF IDF LSI adapted from flyATAC
  
  # binarize the matrix 
  if(binarize){
    message(paste0("Binarizing matrix..."))
    mat@x[mat@x > 0] <- 1 
  }
  
  # order the matrix according to the highly represented genes
  # calculate the number of counts per gene
  if(!is.null(nFeatures)){
    message(paste0("Getting top ", nFeatures, " features..."))
    idx <- head(order(Matrix::rowSums(mat), decreasing = TRUE), nFeatures)
    mat <- mat[idx,] 
  }else{
    idx <- which(Matrix::rowSums(mat) > 0)
    mat <- mat[idx,]
  }
  
  #Calc RowSums and ColSums
  colSm <- Matrix::colSums(mat) # library size per cell
  rowSm <- Matrix::rowSums(mat) # counts per gene 
  
  #Calc TF IDF
  # Term frequency and inverse document frequency
  message("Computing Term Frequency IDF...")
  freqs <- t(t(mat)/colSm) # calculate the term frequency 
  idf   <- as(log(1 + ncol(mat) / rowSm), "sparseVector") # inverse document frequency
  tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs # matrix multiplication %*%
  
  #Calc SVD then LSI 
  # SVD = Singular value decomposition (alternative to PCA)
  message("Computing SVD using irlba...")
  svd <- irlba::irlba(tfidf, nComponents, nComponents)
  svdDiag <- matrix(0, nrow=nComponents, ncol=nComponents)
  diag(svdDiag) <- svd$d
  matSVD <- t(svdDiag %*% t(svd$v))
  rownames(matSVD) <- colnames(mat)
  colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))
  
  #Return Object with all variables
  out <- list(
    matSVD = matSVD, 
    rowSm = rowSm, 
    colSm = colSm, 
    idx = idx, 
    svd = svd, 
    binarize = binarize, 
    nComponents = nComponents,
    date = Sys.Date(),
    seed = 1)
  
  out
  
}

#' **Find nearest neighbors with help of seurat**
#' Use this function to perform PCA on the SVD data matrix and
#' further find the nearest neigbors based on that. From this, communities
#' of cells can be defined as cluster
#' seuratSNN(
#'      matSVD = matSVD (input is a dimensionality reduced matrix) 
#'      dims.use = 1:50 (these are the dimensions used for PCA)
#' )
#' 
#' 
#Clustering function using seurat SNN
seuratSNN <- function(matSVD, dims.use = 1:50, ...){
  set.seed(1)
  message("Making Seurat Object...")
  
  # get original order of rows
  rn <- rownames(matSVD)
  
  # subset the matrix to three
  tmp <- matrix(rnorm(nrow(matSVD) * 3, 10), ncol = nrow(matSVD), nrow = 3)
  colnames(tmp) <- rownames(matSVD)
  rownames(tmp) <- paste0("t",seq_len(nrow(tmp)))
  
  # create seurat object and runn PCA on the SVD matrix
  obj <- Seurat::CreateSeuratObject(tmp, project='scATAC', min.cells=0, min.features=0)
  obj[['pca']] <- Seurat::CreateDimReducObject(embeddings=matSVD, key='PC_', assay='RNA')
  clustParams$object <- obj
  clustParams$reduction <- "pca"
  clustParams$dims <- seq_len(ncol(matSVD))
  
  obj <- suppressWarnings(do.call(Seurat::FindNeighbors, clustParams))
  clustParams$object <- obj
  
  obj <- suppressWarnings(do.call(Seurat::FindClusters, clustParams))
  
  #Get Output
  clust <- obj@meta.data[,ncol(obj@meta.data)]
  clust <- paste0("Cluster",match(clust, unique(clust)))
  names(clust) <- rownames(matSVD)
  clust <- clust[rn]
  
  clust
  
  
}

#' **Generate a sparse matrix from a fragment file**
#' Use this function to create a GRange object that contains chr, start, end, as well as 
#' information about the GC and AT content of each window.
#' countInsertions(
#'      windows = windows (GRange object resulting from makeWindows function), 
#'      fragments = , 
#'      by = "RG"
#' )
#' 

#Sparse Variances Rcpp
sourceCpp(code='
  #include <Rcpp.h>

  using namespace Rcpp;
  using namespace std;

  // [[Rcpp::export]]
  Rcpp::NumericVector computeSparseRowVariances(IntegerVector j, NumericVector val, NumericVector rm, int n) {
    const int nv = j.size();
    const int nm = rm.size();
    Rcpp::NumericVector rv(nm);
    Rcpp::NumericVector rit(nm);
    int current;
    // Calculate RowVars Initial
    for (int i = 0; i < nv; ++i) {
      current = j(i) - 1;
      rv(current) = rv(current) + (val(i) - rm(current)) * (val(i) - rm(current));
      rit(current) = rit(current) + 1;
    }
    // Calculate Remainder Variance
    for (int i = 0; i < nm; ++i) {
      rv(i) = rv(i) + (n - rit(i))*rm(i)*rm(i);
    }
    rv = rv / (n - 1);
    return(rv);
  }'
)

#Compute Fast Sparse Row Variances
sparseRowVariances <- function (m){
  rM <- Matrix::rowMeans(m)
  rV <- computeSparseRowVariances(m@i + 1, m@x, rM, ncol(m))
  return(rV)
}

#' **Generate a sparse matrix from a fragment file**
#' Use this function to create a GRange object that contains chr, start, end, as well as 
#' information about the GC and AT content of each window.
#' countInsertions(
#'      windows = windows (GRange object resulting from makeWindows function), 
#'      fragments = , 
#'      by = "RG"
#' )
#' 

#Helper function for summing sparse matrix groups
groupSums <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups) == ncol(mat))
  gm <- lapply(unique(groups), function(x) {
    if (sparse) {
      Matrix::rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
    else {
      rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
  }) %>% Reduce("cbind", .)
  colnames(gm) <- unique(groups)
  return(gm)
}

#' **Generate a sparse matrix from a fragment file**
#' Use this function to create a GRange object that contains chr, start, end, as well as 
#' information about the GC and AT content of each window.
#' countInsertions(
#'      windows = windows (GRange object resulting from makeWindows function), 
#'      fragments = , 
#'      by = "RG"
#' )
#' 

#Optimized LSI for scRNA-seq analysis
optimizeLSI <- function(mat, scaleTo = 10000, priorCount = 3, pcsUse = 1:25, 
                        resolution = c(0.2, 0.4, 0.8), varFeatures = c(2500, 2500, 2500), seed = 1){
  
  set.seed(seed)
  stopifnot(length(resolution) > 1)
  
  #Initialize List
  lsiOut <- list()
  
  #Initial LSI uses variances that are across all single cells and will have larger batch relationships
  i <- 1
  message("Initial LSI...")
  matNorm <- t(t(mat)/Matrix::colSums(mat)) * scaleTo
  matNorm@x <- log2(matNorm@x + 1)
  idVarFeatures <- head(order(sparseRowVariances(matNorm),decreasing=TRUE), varFeatures[i])
  lsiObj <- calcLSI(mat[idVarFeatures,], binarize = FALSE, nComponents = max(pcsUse))
  clusters <- seuratSNN(lsiObj$matSVD, dims.use = pcsUse, resolution = resolution[i], n.start = 10, print.output = FALSE)
  
  #Store
  lsiOut[[paste0("iter", i)]] <- list(
    lsiMat = lsiObj$matSVD, 
    varFeatures = idVarFeatures, 
    clusters = clusters
  )
  
  for(i in seq(2, length(varFeatures))){
    
    message(sprintf("Additional LSI %s...", i))
    
    #Run LSI
    clusterMat <- edgeR::cpm(groupSums(mat, clusters, sparse = TRUE), log=TRUE, prior.count = priorCount)
    idVarFeatures <- head(order(rowVars(clusterMat), decreasing=TRUE), varFeatures[i])
    lsiObj <- calcLSI(mat[idVarFeatures,], binarize = FALSE, nComponents = max(pcsUse))
    clusters <- seuratSNN(lsiObj$matSVD, dims.use = pcsUse, resolution = resolution[i], n.start = 10, print.output = FALSE)
    
    if(i == length(varFeatures)){
      #Save All Information from LSI Attempt
      lsiOut[[paste0("iter", i)]] <- list(
        lsiObj = lsiObj, 
        varFeatures = idVarFeatures, 
        clusters = clusters,
        matNorm = matNorm
      )
    }else{
      lsiOut[[paste0("iter", i)]] <- list(
        lsiMat = lsiObj$matSVD, 
        varFeatures = idVarFeatures, 
        clusters = clusters
      )
    }
    
  }
  
  return(lsiOut)
  
}

`%notin%` <- Negate(`%in%`)
######################################################################################
# SET WORKING DIRECTORY AND LOAD BATCH CORRECTED (COMBAT) ORGANOID DATA
######################################################################################
# working directory
w.dir <- "/labs/ccurtis2/mjprzy/scRNA_analysis/hashECB_data_freeze/"

# output directory
o.dir <- "/labs/ccurtis2/mjprzy/scRNA_analysis/LSI_projection/"
dir.create(o.dir)
setwd(o.dir)

# read in the COMBAT corrected organoid data
EL.bc.matrix <- read.table("/labs/ccurtis2/mjprzy/scRNA_analysis/combat/freeze/EML_Seq_freeze/EML_Seq_freeze_organoid_bc_data.txt", header = T, sep = "\t") # EL data

######################################################################################
# LOAD THE SATHE REFERENCE DATASET THAT IS USED FOR THE PROJECTION
######################################################################################

## SATHE DATASET
# set sample name
sample.tmp <- "merged_data_SatheEtal_EML"
dir.create(sample.tmp)

# read in the seurat object for the annotated Sathe dataset
# this only comprises epithelial cells
seurat.obj.big <- readRDS(paste0(w.dir, "merged_data_SatheEtal/merged_data_SatheEtal_Epithelial_CellType_seurat_obj.rds"))
mean(seurat.obj.big$nFeature_RNA)
reference.data <- seurat.obj.big@assays$RNA@counts

# HERE I REMOVE THE CELLS WHICH TURNED OUT TO BE NON-EPITHELIAL CELLS IN THE CLUSTERING (SEE CODE BELOW)
# remove cells which are not epithelial cells
remove.cells <- read.table("/labs/ccurtis2/mjprzy/scRNA_analysis/LSI_projection/Sathe_cells_to_remove.txt", header = F, col.names = "Cell_barcodes", sep = "\t")
reference.data <- reference.data[, colnames(reference.data) %notin% remove.cells$Cell_barcodes]

# subset the seurat object as well
valid.cells <- colnames(reference.data)
seurat.obj.big <- subset(seurat.obj.big, cells = valid.cells)

######################################################################################
# CONVERT THE REFERENCE DATASET OF INTEREST INTO A SUMMARIZED EXPERIMENT OBJECT
######################################################################################

# get the reference metadata
seurat.metadata <- seurat.obj.big@meta.data

# intersect both matrices to genes in common
gU <- intersect(rownames(reference.data), rownames(EL.bc.matrix))
gU <- gU[!grepl("^MT", gU)]

# subset the reference dataset to the gene universe
reference.data <- reference.data[gU,]

# create the Summarized Experiment from the reference data for the projection
sce <- SingleCellExperiment(list(counts=reference.data),
                            colData=DataFrame(Cell_barcode = colnames(reference.data)),
                            rowData=DataFrame(hgnc_symbol = rownames(reference.data)),
                            metadata=seurat.metadata
)

# check object
sce

## SATHE DATASET
# class: SingleCellExperiment 
# dim: 15514 10598 
# metadata(21): orig.ident nCount_RNA ... seurat_clusters CellType_epithelial
# assays(1): counts
# rownames(15514): LINC00115 FAM41C ... LINC01669 TSPEAR
# rowData names(1): hgnc_symbol
# colnames(10598): 5846_n1_AAACCTGCAATGGATA-1 5846_n1_AAAGATGGTTCCACAA-1 ... 6709_t1_TCTTTCCCATCACGAT-1 6709_t1_TTCTACACAATGACCT-1
# colData names(1): Cell_barcode
# reducedDimNames(0):
#   spikeNames(0):
#   altExpNames(0):

# class: SingleCellExperiment 
# dim: 15482 6080 
# metadata(23): orig.ident nCount_RNA ... SCT_snn_res.0.5 Celltype
# assays(1): counts
# rownames(15482): LINC00115 FAM41C ... HBZ EMILIN3
# rowData names(1): hgnc_symbol
# colnames(6080): 5846_n1_AAATGCCTCTATCCCG-1 5846_n1_AACTCCCAGCATCATC-1 ... 6709_t1_TCTTTCCCATCACGAT-1 6709_t1_TTCTACACAATGACCT-1
# colData names(1): Cell_barcode
# reducedDimNames(0):
#   spikeNames(0):
#   altExpNames(0):

######################################################################################
# PERFORM ITERATIVE LSI ON THE SCRNA-SEQ DATA
######################################################################################

# initialize the parameters
nPCs <- 1:25 #Number of PCs for clustering
nTop <- c(3000, 3000, 3000) #Choose a higher number of variable peaks
resolution <- c(0.2,0.5,0.8) #Clustering resolutions for Seurat SNN
verbose = TRUE
tstart <- NULL

# set clustering parameters for seurats find cluster method
clustParams <- list()
clustParams$verbose <- verbose
clustParams$tstart <- tstart

# run iterativeLSI procedure on the scRNA data
lsiObj <- optimizeLSI(assay(sce), 
                      resolution = resolution, 
                      pcsUse = nPCs,
                      varFeatures = nTop)

# add the iterativeLSI results
metadata(sce)$optimizeLSI <- lsiObj

# save the SVD matrix
metadata(sce)$matSVD <- lsiObj[[length(lsiObj)]][[1]][[1]] #Last one

# the variable genes
metadata(sce)$variableGenes <- rownames(sce)[lsiObj[[length(lsiObj)]]$varFeatures]

# add iterativeLSI clusters to the 
colData(sce)$Clusters <- lsiObj$iter3$clusters

######################################################################################
# CREATE UMAP EMBEDDING BASED ON ITERATIVE LSI 3
######################################################################################
# set the SVD matrix and the respective clusters from iterativeLSI
matSVD <- metadata(sce)$matSVD
clusters <- colData(sce)$Clusters

# Set Seed and perform UMAP on SVD Matrix
set.seed(1)
uwotUmap <- uwot::umap(
  matSVD, 
  n_neighbors = 35, 
  min_dist = 0.45, 
  metric = "euclidean", 
  n_threads = 1, 
  verbose = TRUE, 
  ret_nn = TRUE,
  ret_model = TRUE
)

# plot the umap embedding
pdf(paste0(sample.tmp, "/", sample.tmp, "_Plot_UMAP-NN-35-MD-45.pdf"), width = 12, height = 12, useDingbats = FALSE)
df <- data.frame(
  x = uwotUmap[[1]][,1],
  y = uwotUmap[[1]][,2], 
  color = clusters
)
ggplot(df,aes(x,y,color=color)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values=viridis(length(unique(clusters)))) +
  xlab("UMAP Dimension 1") + 
  ylab("UMAP Dimension 2")
dev.off()

# Add UMAP coordinates to column data in summarized experiment
colData(sce)$UMAP1 <- uwotUmap[[1]][,1]
colData(sce)$UMAP2 <- uwotUmap[[1]][,2]

# Add UMAP Params for documentation
metadata(sce)$UMAP_Params <- list(NN = 35, MD = 0.45, PCs = 1:25, VarGenes = 3000, Res = "2.6.1")

# Save Summarized Experiment
saveRDS(sce, paste0(sample.tmp, "/SummarizedExperiment_", sample.tmp, ".rds"))

# Save UMAP embedding
save_uwot(uwotUmap, paste0(o.dir, "/", sample.tmp, "/SummarizedExperiment_",sample.tmp, "_UMAP-model.uwot"))

######################################################################################
# PERFORM CELL TYPE ASSIGNMENT FOR THE ITERATIVE LSI CLUSTERS
######################################################################################

# add the normalized count matrix to the seurat object
seurat.obj.big[["LSI"]] <- CreateAssayObject(data = lsiObj$iter3$matNorm)

# switch default assay to LSI
DefaultAssay(object = seurat.obj.big) <- "LSI"

# add the pca and iterative LSI coordinates here
seurat.obj.big[["LSI_pca"]] <- CreateDimReducObject(embeddings = lsiObj$iter3$lsiObj$matSVD, key = "LSIpca_", assay = "LSI")
seurat.obj.big[["LSI_umap"]] <- CreateDimReducObject(embeddings = as.matrix(colData(sce)[,c(3,4)]), key = "LSIumap_", assay = "LSI")

# Find variable features for the RTAM assay
seurat.obj.big <- FindVariableFeatures(object = seurat.obj.big, assay = "LSI", nfeatures = 3000)

# add the clusters to the metadata
seurat.obj.big <- AddMetaData(seurat.obj.big, lsiObj$iter3$clusters, col.name = "iterativeLSI_clusters")
seurat.obj.big$iterativeLSI_clusters <- str_split_fixed(seurat.obj.big$iterativeLSI_clusters, "Cluster", 2)[,2]
seurat.obj.big$iterativeLSI_clusters <- factor(seurat.obj.big$iterativeLSI_clusters, levels = c(1:length(unique(seurat.obj.big$iterativeLSI_clusters))))

# check the clustering representation
pdf(paste0(sample.tmp, "/", sample.tmp, "_CellType_LSI_umap_plot.pdf"), width = 20, height = 20)
print(DimPlot(seurat.obj.big, reduction = "LSI_umap", label = T, pt.size = 1))
print(DimPlot(seurat.obj.big, reduction = "LSI_umap", group.by = "orig.ident", pt.size = 1))
dev.off()

pdf(paste0(sample.tmp, "/", sample.tmp, "_EpithelialMarkers_LSI_umap_Featureplot.pdf"), width = 20, height = 20)
FeaturePlot(seurat.obj.big, reduction = "LSI_umap", features = c("EPCAM", "KRT18", "MUC1", "KRT19", "CDH1", "CLDN4"))
FeaturePlot(seurat.obj.big, reduction = "LSI_umap", features = c("CD4", "VIM", "ACTA2", "PTPRC"))
dev.off()

pdf(paste0(sample.tmp, "/", sample.tmp, "_IntestinalGenes_Sathe_LSI_umap_Featureplot.pdf"), width = 20, height = 20)
FeaturePlot(seurat.obj.big, reduction = "LSI_umap", features = c("TFF3", "FABP1", "SPINK4", "MUC13", "REG4", "SOX4", "HES1"))
dev.off()

# change cell identities for new assignment based on iterative LSI#
seurat.obj.big$epithelialOnly_CellType <- Idents(seurat.obj.big)
Idents(object = seurat.obj.big) <- "iterativeLSI_clusters"

# plot probability distributions across clusters for the top5 genes indicated in the Supplementary data from Zhang et al.
# However, markers shared by different cell types were removed from both.
pdf(paste0(sample.tmp, "/", sample.tmp, "_exp_violin_plot.pdf"), width = 20, height = 20)
# check for PMCs
print(VlnPlot(object = seurat.obj.big, assay = "LSI", features = c("GKN1", "GKN2", "MUC5AC", "TFF1", "DPCR1")) + labs(title = "PMCs"))
# check for MSCs
print(VlnPlot(object = seurat.obj.big, assay = "LSI", features = c("OLFM4", "RPL7", "CLDN4","TSPAN8", "REG1A")) + labs(title = "MSCs"))
# check for enterocytes
print(VlnPlot(object = seurat.obj.big, assay = "LSI", features = c("FABP1", "FABP2", "RBP2", "ANPEP", "APOA4")) + labs(title = "Enterocytes"))
# check for GMC
print(VlnPlot(object = seurat.obj.big, assay = "LSI", features = c("PRR4", "C6orf58","MUC6", "TFF2", "LTF")) + labs(title = "GMCs"))
# check for enteroendocrine
print(VlnPlot(object = seurat.obj.big, assay = "LSI", features = c("CHGA", "PCSK1N", "SCG5", "CHGB", "TPH1")) + labs(title = "Enteroendocrine")) 
# check for PCs
print(VlnPlot(object = seurat.obj.big, assay = "LSI", features = c("TOP2A", "MKI67", "UBE2C", "HMGB2", "PTTG1")) + labs(title = "PCs"))
# check for cancer cells
print(VlnPlot(object = seurat.obj.big, assay = "LSI", features = c("REG4","CLDN7","KRT18", "LGALS3", "CEACAM6")) + labs(title = "Cancer Cells"))
print(VlnPlot(object = seurat.obj.big, assay = "LSI", features = c("CLDN3", "CST1", "MUC3A", "CLDN4", "PI3", "UBD")) + labs(title = "Cancer Cells"))
print(VlnPlot(object = seurat.obj.big, assay = "LSI", features = c("CDH17", "PRAP1", "UBE2C", "CCL20", "LCN2", "SERPINB5")) + labs(title = "Cancer Cells"))
print(VlnPlot(object = seurat.obj.big, assay = "LSI",  features = c("RRM2", "MYBL2", "MMP7", "TPX2", "MISP", "TMPRSS4")) + labs(title = "Cancer Cells"))
print(VlnPlot(object = seurat.obj.big, assay = "LSI", features = c("RMRP", "CLDN1", "GPRC5A", "CLRN3", "CXCL1", "MSLN")) + labs(title = "Cancer Cells"))
# check for neck like cells
print(VlnPlot(object = seurat.obj.big, assay = "LSI", features = c("MIA", "CXCL3", "CXCL2", "CXCL17", "CLU")) + labs(title = "Neck-like cells"))
# check for Goblet cells
print(VlnPlot(object = seurat.obj.big, assay = "LSI", features = c("SPINK4", "TFF3", "MUC2", "ITLN1", "ZG16")) + labs(title = "Goblet Cells"))
# check for Chief cells
print(VlnPlot(object = seurat.obj.big, assay = "LSI", features = c("PGA3", "PGA4", "LIPF", "CHIA", "PDIA2")) + labs(title = "Chief Cells")) 
# check for parietal cells
print(VlnPlot(object = seurat.obj.big, assay = "LSI", features = c("ATP4A", "ATP4B", "VEGFB")) + labs(title = "Parietal cells"))
# check for endocrine cells
print(VlnPlot(object = seurat.obj.big, assay = "LSI", features = c("GAST", "GHRL", "SST")) + labs(title = "Endocrine cells"))
# check for Metaplasia-like Cells
print(VlnPlot(object = seurat.obj.big, assay = "LSI", features = c("REG4", "ACE2", "PCK1", "RBP2", "MTTP", "SLC26A3")) + labs(title = "Metaplasia-like Cells"))
print(VlnPlot(object = seurat.obj.big, assay = "LSI", features = c("SI", "ANPEP", "APOB", "CPS1", "GKN2", "S100P")) + labs(title = "Metaplasia-like Cells"))
print(VlnPlot(object = seurat.obj.big, assay = "LSI", features = c("FCGBP", "LGALS4", "GDA", "LYZ", "CFTR", "KRT20")) + labs(title = "Metaplasia-like Cells"))
print(VlnPlot(object = seurat.obj.big, assay = "LSI", features = c("ADH1C", "AKR1B10", "CDCA7", "SLC5A1", "CYP2C18", "ELOVL6")) + labs(title = "Metaplasia-like Cells"))
print(VlnPlot(object = seurat.obj.big, assay = "LSI", features = c("MUC13", "SLC6A14", "AADAC", "HSD17B2", "GCNT3")) + labs(title = "Metaplasia-like Cells"))
dev.off()

# get all cluster markers to compare clusters individually 
seurat.obj.big.markers <- FindAllMarkers(seurat.obj.big, assay = "LSI", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
seurat.obj.big.markers <- seurat.obj.big.markers[order(seurat.obj.big.markers$cluster, seurat.obj.big.markers$avg_logFC,decreasing = c(F,T)),]
write.table(seurat.obj.big.markers, paste0(sample.tmp, "/", sample.tmp, "_clusterMarkers_0.25pct_0.5logFC.txt"), sep = "\t", quote = F, col.names = T)

# HERE WE SEE THAT SOME CELLS WERE NON-EPITHLIAL WHICH IS WHY THEY WERE EXCLUDED IN A SECOND ITERATION
# cells.to.remove <- rownames(seurat.obj.big@meta.data[seurat.obj.big$iterativeLSI_clusters == 8,])
# write.table(cells.to.remove, "/labs/ccurtis2/mjprzy/scRNA_analysis/LSI_projection/Sathe_cells_to_remove.txt", quote = F, col.names = F, row.names = F)

######################################################################################
# CHECK AND REASSIGN CLUSTERS FOR SATHE ET AL
######################################################################################
## SATHE DATASET
# ASSIGN CELL TYPES TO THE CLUSTERS FOR SATHE DATASET
new.cluster.ids <- c("Chief Cells", "PMCs", "Chief Cells", "Enterocytes", "PMCs", "Parietal Cells", "Enteroendocrine Cells", "Enteroendocrine Cells", 
                     "Goblet Cells", "Goblet Cells", "PMCs", "MSCs", "GMCs", "Mucosal-like malignant Cells",
                     "Enterocytes", "Mucosal-like malignant Cells", "Non-Mucosal-like malignant Cells", "Proliferative Cells", "Non-Mucosal-like malignant Cells", "Proliferative Cells",
                     "Enterocytes", "Paneth Cells")

# change cluster ids for Cell names
names(new.cluster.ids) <- levels(seurat.obj.big)
seurat.obj.big <- RenameIdents(seurat.obj.big, new.cluster.ids)

# save to metadata as well
seurat.obj.big$LSI_CellType <- Idents(seurat.obj.big)

######################################################################################
# ANALYSE THE DIFFERENCE BETWEEN THE NON-MUCOSAL & MUCOSAL-LIKE MALIGNANT CELL CLUSTER
######################################################################################

# compare cluster to all other malignant cells
malignant.markers <- FindMarkers(seurat.obj.big, assay = "LSI", ident.1 = "Non-Mucosal-like malignant Cells", ident.2 = "Mucosal-like malignant Cells", min.pct = 0.5, logfc.threshold = 0.5)
malignant.markers <- malignant.markers[order(malignant.markers$avg_logFC, decreasing = T),]
write.table(malignant.markers, paste0(sample.tmp, "/", sample.tmp, "_MalignantClusters_0.5pct_0.5logFC_DESeq2.txt"), sep = "\t", quote = F, col.names = T)

######################################################################################
# QUANTIFY THE ABUNDANCE OF SAMPLES PER CLUSTER/CELL TYPE
######################################################################################

## SAMPLE ID VS CELL TYPE
# get the contribution of each sample to each cell type
contingency.table <- table(seurat.obj.big$orig.ident, seurat.obj.big$LSI_CellType)
contingency.table <- round(prop.table(contingency.table, 2), 3)
contingency.table <- as.data.frame.matrix(contingency.table)
contingency.table <- as.matrix(contingency.table)

# generate the second heatmap with rownames annotated in contrast to the first
ht <- Heatmap(contingency.table,  
              col = circlize::colorRamp2(c(0, 1), c("white", "darkblue")), 
              cluster_rows = FALSE, 
              cluster_columns = FALSE, 
              name = "High_threshold_pred", 
              row_names_gp = gpar(fontsize = 15, fontface = "bold"),
              column_names_gp = gpar(fontsize = 15, fontface = "bold"),
              show_column_names = TRUE,
              column_names_side = "bottom",
              column_names_rot = 45,
              show_row_names = TRUE,
              rect_gp = gpar(col = "black", lwd = 1), 
              column_title_gp = gpar(fontsize = 15, fontface = "bold"), 
              row_title_gp = gpar(fontsize = 15, fontface = "bold"), 
              border = T,
              cell_fun = function(j, i, x, y, width, height, fill) { # j: column index in matrix, i: row index, x: x coordinate of cell middle point, y: y coordinate for cell middpoint, 
                grid.text(sprintf("%.2f", contingency.table[i, j]), x, y, gp = gpar(fontsize = 10))
              },
              column_title = "Cell types", 
              row_title = "Sample ID")

# plot and save the list as a complete figure
pdf(paste0(sample.tmp, "/", sample.tmp, "_ContingencyTable_SampleID_CellType.pdf"), width = 15, height = 15)
print(ht)
dev.off()

## SAMPLE ID VS CLUSTER
# get the contribution of each sample to each cluster
contingency.table <- table(seurat.obj.big$iterativeLSI_clusters, seurat.obj.big$orig.ident)
contingency.table <- round(prop.table(contingency.table, 2), 3)
contingency.table <- as.data.frame.matrix(contingency.table)
contingency.table <- as.matrix(contingency.table)

# generate the second heatmap with rownames annotated in contrast to the first
ht <- Heatmap(contingency.table,  
              col = circlize::colorRamp2(c(0, 1), c("white", "darkblue")), 
              cluster_rows = FALSE, 
              cluster_columns = FALSE, 
              name = "High_threshold_pred", 
              row_names_gp = gpar(fontsize = 15, fontface = "bold"),
              column_names_gp = gpar(fontsize = 15, fontface = "bold"),
              show_column_names = TRUE,
              column_names_side = "bottom",
              column_names_rot = 45,
              show_row_names = TRUE,
              rect_gp = gpar(col = "black", lwd = 1), 
              column_title_gp = gpar(fontsize = 15, fontface = "bold"), 
              row_title_gp = gpar(fontsize = 15, fontface = "bold"), 
              border = T,
              cell_fun = function(j, i, x, y, width, height, fill) { # j: column index in matrix, i: row index, x: x coordinate of cell middle point, y: y coordinate for cell middpoint, 
                grid.text(sprintf("%.2f", contingency.table[i, j]), x, y, gp = gpar(fontsize = 10))
              },
              column_title = "Sample ID", 
              row_title = "Clusters")

# plot and save the list as a complete figure
pdf(paste0(sample.tmp, "/", sample.tmp, "_ContingencyTable_SampleID_Cluster.pdf"), width = 15, height = 15)
print(ht)
dev.off()

######################################################################################
# PLOT UMAP EMBEDDINGS FOR THE NEW REPRESENTATION
######################################################################################

## SAMPLE TYPE (HISTOLOGICAL IDENTITY)
plot1 <- DimPlot(seurat.obj.big, reduction = "LSI_umap", group.by = "sample_type", pt.size = 1.5) +
  xlab("UMAP Dimension 1") + 
  ylab("UMAP Dimension 2") +
  theme_classic() +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=25, face = "bold")) + 
  theme(legend.text = element_text(colour="black", size=25, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  guides(color=guide_legend("Sample types", override.aes = list(size = 6)))

# save single plot
pdf(paste0(sample.tmp, "/", sample.tmp, "_iterativeLSI_UMAP_sampleType.pdf"), width = 6.5, height = 6)
print(plot1 + NoLegend() + NoAxes()) 
dev.off()

# Using the cowplot package
legend <- cowplot::get_legend(plot1)

grid.newpage()
grid.draw(legend)

## ITERATIVE LSI CLUSTERS
seurat.obj.big$iterativeLSI_clusters <- factor(seurat.obj.big$iterativeLSI_clusters, levels = c(1:length(unique(seurat.obj.big$iterativeLSI_clusters))))
plot2 <- DimPlot(seurat.obj.big, reduction = "LSI_umap", group.by = "iterativeLSI_clusters", label = T, label.size = 5, pt.size = 1.5) +
  xlab("UMAP Dimension 1") + 
  ylab("UMAP Dimension 2") +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=25, face = "bold")) + 
  theme(legend.text = element_text(colour="black", size=25, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  guides(color=guide_legend("Clusters", override.aes = list(size = 8)))

# save single plot
pdf(paste0(sample.tmp, "/", sample.tmp, "_iterativeLSI_UMAP_clusters.pdf"), width = 6.5, height = 6)
print(plot2 + NoLegend() + NoAxes()) 
dev.off()

# Using the cowplot package
legend <- cowplot::get_legend(plot2)

grid.newpage()
grid.draw(legend)

## ANNOTATED CELL TYPE IDENTITY
plot3 <- DimPlot(seurat.obj.big, reduction = "LSI_umap", label = T, label.size = 3, pt.size = 1.5) +
  xlab("UMAP Dimension 1") + 
  ylab("UMAP Dimension 2") +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=25, face = "bold")) + 
  theme(legend.text = element_text(colour="black", size=25, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  guides(color=guide_legend("Cell Types", override.aes = list(size = 6)))

# save single plot
pdf(paste0(sample.tmp, "/", sample.tmp, "_iterativeLSI_UMAP_CellTypes.pdf"), width = 6.5, height = 6)
print(plot3 + NoLegend() + NoAxes()) 
dev.off()

# Using the cowplot package
legend <- cowplot::get_legend(plot3)

grid.newpage()
grid.draw(legend)

# make plotlist
plotlist <- list(plot1, plot2, plot3)

# plot the list
nCol <- 3
pdf(paste0(sample.tmp, "/", sample.tmp, "_iterativeLSI_Manual_Comparison_plot.pdf"), width = 60, height = 20)
cowplot::plot_grid(plotlist = plotlist, ncol = nCol)
dev.off()

# switch default back to SCT
DefaultAssay(object = seurat.obj.big) <- "SCT"

# save the seurat object with the iterativeLSI assay
saveRDS(seurat.obj.big, paste0(sample.tmp, "/", sample.tmp, "_iterativeLSI_seurat.obj.rds"))
seurat.obj.big <- readRDS(paste0(sample.tmp, "/", sample.tmp, "_iterativeLSI_seurat.obj.rds"))
table(seurat.obj.big$sample_type)
      
