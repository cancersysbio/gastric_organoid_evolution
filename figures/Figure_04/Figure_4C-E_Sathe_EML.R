#############################################################################################################################
##                                                                                                                      
##  Perform clustering and umap embedding on 10x single-cell RNA-sequencing data Clustering and scRNA-seq UMAP for EML organoid data
##
##  Date: 16 May 2021                                                                                                                   
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
list.of.packages <- c("reshape2", "optparse", "BSgenome", "RColorBrewer", "ggplot2", "scales",  "dendextend", "tidyverse", 
                      "Matrix", "ComplexHeatmap", "Rtsne", "robustbase", "psych", "cluster", "Matrix.utils",
                      "BSgenome.Hsapiens.UCSC.hg38", "Seurat", "GenomicRanges", "SingleCellExperiment", "matrixStats", "readr",
                      "magrittr", "edgeR", "uwot", "BiocManager", "Rcpp", "biomaRt", "dplyr", "viridis", "FNN", "httr", "Matrix","reshape2",
                      "ggrepel")
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

projectLSI <- function(mat, lsi){   
  
  #Get Same Features
  mat <- mat[lsi$idx,]
  if(lsi$binarize){
    message(paste0("Binarizing matrix..."))
    mat@x[mat@x > 0] <- 1       
  }
  
  #Calc TF IDF
  rowsToZero <- which(lsi$rowSm == 0)
  setToZero <- which((mat@i + 1) %in% rowsToZero)
  if(length(setToZero) > 0){
    mat@x[setToZero] <- 0
  }
  
  message("Computing Term Frequency IDF...")
  freqs <- t(t(mat)/Matrix::colSums(mat))
  idf   <- as(log(1 + length(lsi$colSm) / lsi$rowSm), "sparseVector")
  tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs
  if(length(Matrix::which(is.na(tfidf),arr.ind=TRUE)) > 0){
    tfidf[Matrix::which(is.na(tfidf),arr.ind=TRUE)] <- 0 #weird Inf * 0
  }
  
  #Calc V
  V <- t(tfidf) %*% lsi$svd$u %*% diag(1/lsi$svd$d)
  
  #Calc SVD then LSI
  message("Computing SVD using irlba...")
  svdDiag <- matrix(0, nrow=lsi$nComponents, ncol=lsi$nComponents)
  diag(svdDiag) <- lsi$svd$d
  matSVD <- t(svdDiag %*% t(V))
  rownames(matSVD) <- colnames(mat)
  colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))
  
  return(matSVD)
  
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

sparseMatTTest <- function(mat1, mat2, m0 = 0){
  #Get Population Values
  n1 <- ncol(mat1)
  n2 <- ncol(mat2)
  n <- n1 + n2
  #Sparse Row Means
  m1 <- Matrix::rowMeans(mat1, na.rm=TRUE)
  m2 <- Matrix::rowMeans(mat2, na.rm=TRUE)
  #Sparse Row Variances
  v1 <- computeSparseRowVariances(mat1@i + 1, mat1@x, m1, n1)
  v2 <- computeSparseRowVariances(mat2@i + 1, mat2@x, m2, n2)
  #Calculate T Statistic
  se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*v1 + (n2-1)*v2)/(n1+n2-2) )
  tstat <- (m1-m2-m0)/se
  #tstat <- sqrt((n1 * n2) / n) / sqrt((n1-1)/(n-2)*v1 + (n2-1)/(n-2)*v2)
  pvalue <- 2*pt(-abs(tstat), n - 2)
  fdr <- p.adjust(pvalue, method = "fdr")
  out <- data.frame(fdr = fdr, pval = pvalue, tstat = tstat, mean1 = m1, mean2 = m2, var1 = v1, var2 = v2, n1 = n1, n2 = n2)
  return(out)
}

`%notin%` <- Negate(`%in%`)

is.integer0 <- function(x){
  is.integer(x) && length(x) == 0L
}

######################################################################################
# SET WORKING DIRECTORY AND LOAD BATCH CORRECTED (COMBAT) ORGANOID DATA
######################################################################################
# working directory
w.dir <- "/labs/ccurtis2/mjprzy/scRNA_analysis/hashECB_data_freeze/"

# output directory
o.dir <- "/labs/ccurtis2/mjprzy/scRNA_analysis/LSI_projection/freeze"
dir.create(o.dir)
setwd(o.dir)

# read in the COMBAT corrected organoid data
EML.bc.matrix <- read.table("/labs/ccurtis2/mjprzy/scRNA_analysis/combat/freeze/EML_Seq_freeze/EML_Seq_freeze_organoid_bc_data.txt", header = T, sep = "\t")
colnames(EML.bc.matrix) <- gsub("X", "", colnames(EML.bc.matrix))
EML.bc.matrix <- as.matrix(EML.bc.matrix)
EML.data <- as(EML.bc.matrix, "sparseMatrix")
save(EML.data, file = "/labs/ccurtis2/mjprzy/scRNA_analysis/combat/freeze/EML_Seq_freeze/EML_Seq_freeze_organoid_bc.RData", compress = T)

# read in Seurat object to get the corresponding metadata
seurat.obj.big <- readRDS(paste0(w.dir, "EML_Seq_freeze/EML_Seq_freeze_Adapted_Patient_seurat_obj.rds")) # EL data
seurat.metadata <- seurat.obj.big@meta.data

######################################################################################
# LOAD THE SATHE REFERENCE DATASET THAT IS USED FOR THE PROJECTION
######################################################################################
## SATHE DATASET
# set sample name
sample.tmp <- "EML_Seq_freeze_Sathe"
dir.create(sample.tmp)

# read in the seurat object for the annotated Sathe dataset which was generated 
# within the scRNA_01 script
seurat.obj.big <- readRDS(paste0("/labs/ccurtis2/mjprzy/scRNA_analysis/LSI_projection/merged_data_SatheEtal_EML/merged_data_SatheEtal_EML_iterativeLSI_seurat.obj.rds"))
reference.data <- seurat.obj.big@assays$RNA@counts

######################################################################################
# CONVERT THE REFERENCE DATASET OF INTEREST INTO A SUMMARIZED EXPERIMENT OBJECT
######################################################################################

#Identify Gene Universe
gU <- intersect(rownames(EML.data), rownames(reference.data))
gU <- gU[!grepl("^MT", gU)]

# subset the EL data to the same gene universe
EML.data <- EML.data[gU,]

# create the Summarized Experiment from the normalized data
sce <- SingleCellExperiment(list(counts=EML.data),
                            colData=DataFrame(Cell_barcode = colnames(EML.data)),
                            rowData=DataFrame(hgnc_symbol = rownames(EML.data)),
                            metadata=seurat.metadata
)

# add the patient name here
colData(sce)$Group <- str_split_fixed(rownames(colData(sce)), "_", 2)[,1]

# check Summarized Experiment object
sce

## EL Organoids
# class: SingleCellExperiment 
# dim: 15514 31772 
# metadata(26): Cell_Barcode orig.ident ... DF.classifications_0.25_0.3_2562 DF.ID
# assays(1): counts
# rownames(15514): AL627309.1 LINC00115 ... CLDN14 LINC00114
# rowData names(1): hgnc_symbol
# colnames(31772): 6077_EL_AAACCCAAGACATCAA 6077_EL_AAACCCAAGACGGAAA ... 0891_EL_TTTGTTGTCACACCCT 0891_EL_TTTGTTGTCTAGCCTC
# colData names(2): Cell_barcode Group
# reducedDimNames(0):
#   spikeNames(0):
#   altExpNames(0):

# class: SingleCellExperiment 
# dim: 15482 26145 
# metadata(25): Cell_Barcode orig.ident ... pANN_0.25_0.3_0 Doublet_Score
# assays(1): counts
# rownames(15482): LINC00115 FAM41C ... CLDN14 LINC00114
# rowData names(1): hgnc_symbol
# colnames(26145): 6077_EL_AAACCCAAGACATCAA 6077_EL_AAACCCAAGACGGAAA ... 0891_EL_TTTGTTGTCTAGCCTC 0891_EL_TTTGTTGTCTCAACGA
# colData names(2): Cell_barcode Group
# reducedDimNames(0):
#   spikeNames(0):
#   altExpNames(0):

######################################################################################
# CLUSTERING ANALYSIS OF THE FULL ORGANOID DATASET STARTS HERE
######################################################################################
# Set Clustering Parameters
varGenesToUse <- c(1000, 1000, 1000) #Choose a higher number of variable peaks
resolution <- c(0.2,0.4,0.6) #Clustering resolutions for Seurat SNN
verbose = TRUE
tstart <- NULL

# set clustering parameters for seurats find cluster method
clustParams <- list()
clustParams$verbose <- verbose
clustParams$tstart <- tstart

#Optimize LSI Features for STP samples
matAll <- assay(sce)
lsiObj <- optimizeLSI(matAll, resolution = resolution, varFeatures = varGenesToUse,  pcsUse = 1:25)

#UMAP
set.seed(1)
umap <- uwot::umap(
  lsiObj[[length(lsiObj)]]$lsiObj$matSVD[,1:25], 
  n_neighbors = 35, 
  min_dist = 0.5, 
  metric = "euclidean", 
  n_threads = 5, 
  verbose = TRUE, 
  ret_model = FALSE
)

######################################################################################
# CLASSIFICATION OF EARLY LATE SAMPLES FROM GASTRIC ORGANOIDS
######################################################################################

# add timepoint metadata
sce@metadata$timepoint <- str_split_fixed(sce@metadata$HashTag, "_", 3)[,3]
sce@metadata$timepoint[grep("WT", sce@metadata$HashTag)] <- "WT"

# calculate the percentage of disease cells per cluster Plot Info
cells <- sce@metadata$timepoint

# split the cells by the clusters 
splitCells <- split(cells,lsiObj[[length(lsiObj)]]$clusters)

# calculate the proporationn of LATE cells in each cluster
df <- data.frame(
  clusters = names(splitCells),
  proportion = unlist(lapply(seq_along(splitCells), function(x) sum(splitCells[[x]]=="LATE") / length(splitCells[[x]])))
)

#Plot UMAP Data Frame
plotDF <- data.frame(umap)
rownames(plotDF) <- c(colnames(sce))

# add metadata to the plot data frame
plotDF$type <- cells
plotDF$clusters <- lsiObj[[length(lsiObj)]]$clusters
plotDF$patient_id <- str_split_fixed(rownames(plotDF), "_", 2)[,1]

# setup a projection plot directory
plotDir <- paste0(sample.tmp, "/classification/")
dir.create(plotDir,recursive=TRUE)

####################################################
# LOAD THE SATHE DATASET FROM STEP 1
####################################################

# Here, we read in the data and uwot object step 1
se <- readRDS("/labs/ccurtis2/mjprzy/scRNA_analysis/LSI_projection/merged_data_SatheEtal_EML/SummarizedExperiment_merged_data_SatheEtal_EML.rds")

#Load Saved UMAP Manifold
umapManifold <- uwot::load_uwot("/labs/ccurtis2/mjprzy/scRNA_analysis/LSI_projection/merged_data_SatheEtal_EML/SummarizedExperiment_merged_data_SatheEtal_EML_UMAP-model.uwot")

####################################################
# PERFORM PROJECTION INTO LSI UMAP
####################################################

# subset the organoid LSI Projection Matrix to the variableGenes from the projection
lsiGenes <- metadata(se)$variableGenes
matProjectLSI <- assay(sce[lsiGenes,])

# LSI Project
lsiReference <- metadata(se)$optimizeLSI[[length(metadata(se)$optimizeLSI)]]$lsiObj

# run lsi projection
lsiProjection <- projectLSI(matProjectLSI, lsiReference)

#UMAP Projection
#Set Seed Prior to umap_transform (see uwot github)
set.seed(1)
umapProjection <- uwot::umap_transform(as.matrix(lsiProjection)[,1:25], umapManifold, verbose = TRUE)

######################################################################################
#                   ANALYSIS OF INDIVIDUAL SAMPLES HERE
######################################################################################
# get the sample ids
sample.ids <- unique(metadata(sce)$HashTag)

# alternatively, use patient ids
patient.ids <- unique(str_split_fixed(metadata(sce)$HashTag, "_", 3)[,1])

###############################################################################
# LOAD THE SATHE SEURAT OBJECT WITH THE CELL TYPE ANNOTATION
###############################################################################

# load the respective seurat object for the reference with the assigned cell types
seurat.obj.big <- readRDS(paste0("/labs/ccurtis2/mjprzy/scRNA_analysis/LSI_projection/merged_data_SatheEtal_EML/merged_data_SatheEtal_EML_iterativeLSI_seurat.obj.rds"))
seurat.obj.big$LSI_CellType <- Idents(seurat.obj.big)
table(seurat.obj.big$orig.ident)

#Load Saved UMAP Manifold
umapManifold <- uwot::load_uwot("/labs/ccurtis2/mjprzy/scRNA_analysis/LSI_projection/merged_data_SatheEtal_EML/SummarizedExperiment_merged_data_SatheEtal_EML_UMAP-model.uwot")

######################################################################################
# PLOT UMAP EMBEDDING OF THE REFERENCE/BACKGROUND DATASET ONLY
######################################################################################

# subset to background onlys
pdf(paste0(plotDir, sample.tmp,"_coloured_background_only.pdf"), width = 12, height = 12, useDingbats = FALSE)
DimPlot(seurat.obj.big, reduction = "LSI_umap", pt.size = 1.5) +
  theme_bw() +
  xlab("UMAP Dimension 1") + 
  ylab("UMAP Dimension 2") +
  scale_color_manual(values=c("Chief Cells" = "lightgrey", "Parietal Cells" = "lightgrey", "Enteroendocrine Cells" = "lightgrey", "PMCs" = magma(9)[5], 
                              "Goblet Cells" = "lightgrey", "MSCs" = "lightgrey", "Mucosal-like malignant Cells" = magma(9)[4], "Enterocytes" = "lightgrey",
                              "GMCs" = "lightgrey", "Endocrine Cells" = "lightgrey", "Proliferative Cells" = "lightgrey", "Non-Mucosal-like malignant Cells" = magma(9)[3],
                              "Paneth Cells" = "lightgrey")) +
  theme_classic() + NoLegend() +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=25, face = "bold")) + 
  theme(legend.text = element_text(colour="black", size=25, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  guides(color=guide_legend("Patients", override.aes = list(size = 6))) + NoLegend() + NoAxes()
dev.off()

###############################################################################
# ITERATIVE OVER ALL SAMPLES IN THE EL ORGANOID DATASET
###############################################################################

# set up nearest neighbor dataframe for all samples
nn.table <- table(seurat.obj.big$LSI_CellType) %>% bind_rows
nn.table <- as.data.frame(nn.table)

# setup a projection plot directory
plotDir <- paste0(sample.tmp, "/classification/")
dir.create(plotDir,recursive=TRUE)

i <- 1
# iterative over each individual sample
for (i in 1:length(sample.ids)){
  
  # set id from samples
  id <- sample.ids[i]
  
  # set id from patients
  # id <- patient.ids[i]
  print(id)
  
  # get only the cells associated with this sample
  seDisease <- sce[,grep(id, sce@metadata$HashTag)]
  
  #Identify Gene Universe
  gU <- intersect(rownames(seDisease), rownames(se))
  gU <- gU[!grepl("^MT", gU)]
  
  ######################################################################################
  # CLUSTERING ANALYSIS STARTS HERE
  ######################################################################################
  # Set Clustering Parameters
  varGenesToUse <- c(1000, 1000, 1000) #Choose a higher number of variable peaks
  resolution <- c(0.2,0.8,0.8) #Clustering resolutions for Seurat SNN
  verbose = TRUE
  tstart <- NULL
  input_knn <- 25
  scaleTo <- 10000
  nMax <- 3000
  
  # set clustering parameters for seurats find cluster method
  clustParams <- list()
  clustParams$verbose <- verbose
  clustParams$tstart <- tstart
  
  #Optimize LSI Features for STP samples
  matAll <- cbind(assay(seDisease[gU,]), assay(se[gU,])) # combine single sample with the Zhang et al dataset
  lsiObj <- optimizeLSI(matAll, resolution = resolution, varFeatures = varGenesToUse,  pcsUse = 1:25)
  
  #UMAP
  set.seed(1)
  umap <- uwot::umap(
    lsiObj[[length(lsiObj)]]$lsiObj$matSVD[,1:25], 
    n_neighbors = 35, 
    min_dist = 0.5, 
    metric = "euclidean", 
    n_threads = 5, 
    verbose = TRUE, 
    ret_model = FALSE
  )
  
  ####################################################
  # PERFORM PROJECTION INTO LSI UMAP
  ####################################################
  
  ## in the next step, I guess it is assumed, that the matrices contain mainly overlapping genes
  ## however, this does not seem to be the case, which is why I adapted the objects in the following
  #LSI Projection Matrix
  lsiGenes <- metadata(se)$variableGenes
  matProjectLSI <- assay(seDisease[lsiGenes,])
  
  #LSI Project
  lsiReference <- metadata(se)$optimizeLSI[[length(metadata(se)$optimizeLSI)]]$lsiObj
  
  # run projectionn
  lsiProjection <- projectLSI(matProjectLSI, lsiReference)
  
  #UMAP Projection
  #Set Seed Prior to umap_transform (see uwot github)
  set.seed(1)
  umapProjection <- uwot::umap_transform(as.matrix(lsiProjection)[,1:25], umapManifold, verbose = TRUE)
  
  ######################################################################################
  # PLOT PROJECTION ACCORDING TO SAMPLES
  ######################################################################################
  
  #Plot Projection
  refDF <- data.frame(row.names = colnames(se), X1 = umapManifold$embedding[,1], X2 = umapManifold$embedding[,2], Type = "reference")
  proDF <- data.frame(row.names = colnames(seDisease), X1 = umapProjection[,1], X2 = umapProjection[,2], Type = "organoid")
  projectionDF <- rbind(refDF, proDF)
  
  # get the celltype assignment here + barcode
  celltype.df <- data.frame(Cell_barcodes = rownames(seurat.obj.big@meta.data), CellType = seurat.obj.big$LSI_CellType, clusters = se$Clusters)
  
  # create UMAP dataframe merge with the celltype assignment from the background
  projectionDF$Cell_barcodes <- rownames(projectionDF)
  projectionDF <- merge(projectionDF, celltype.df, by = "Cell_barcodes", all.x = T)
  projectionDF$CellType <- as.character(projectionDF$CellType)
  
  # assign the organoid cell type to the organoids which do not have a celltype assignment
  projectionDF[is.na(projectionDF$CellType), "CellType"] <- "organoid"
  rownames(projectionDF) <- projectionDF$Cell_barcodes
  
  # order the dataframe to project organoids on reference and not the other way round
  projectionDF$Type <- factor(projectionDF$Type, levels = c("reference", "organoid"))
  projectionDF <- projectionDF[order(projectionDF$Type),]
  
  ######################################################################################
  # PLOT PROJECTION ACCORDING TO SAMPLES WITH COLOURING ACCORDING TO EARLY LATE WT
  ######################################################################################
  
  if (!is.integer0(grep("EARLY", id))){
    
    ## EARLY
    pdf(paste0(plotDir,id,"-Projection-UMAP.pdf"), width = 12, height = 12, useDingbats = FALSE)
    print(ggplot(projectionDF, aes(X1,X2,color=Type)) + 
            geom_point() +
            theme_bw() +
            xlab("UMAP Dimension 1") + 
            ylab("UMAP Dimension 2") +
            scale_color_manual(values=c("reference"="lightgrey", "organoid" = "#EA9F37")) +
            theme_classic() +
            theme(axis.line=element_blank(),axis.text.x=element_blank(),
                  axis.text.y=element_blank(),axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),legend.position="none",
                  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),plot.background=element_blank()) +
            geom_density_2d(data=projectionDF[projectionDF$Type == "organoid",], mapping=aes(x=X1,y=X2),color="#EA9F37",n=300,h=3)) + NoLegend() + NoAxes()
    dev.off()
    
  } else if (!is.integer0(grep("LATE", id))){
    
    ## LATE
    pdf(paste0(plotDir,id,"-Projection-UMAP.pdf"), width = 12, height = 12, useDingbats = FALSE)
    print(ggplot(projectionDF, aes(X1,X2,color=Type)) + 
            geom_point() +
            theme_bw() +
            xlab("UMAP Dimension 1") + 
            ylab("UMAP Dimension 2") +
            scale_color_manual(values=c("reference"="lightgrey", "organoid" = "#782867")) +
            theme_classic() +
            theme(axis.line=element_blank(),axis.text.x=element_blank(),
                  axis.text.y=element_blank(),axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),legend.position="none",
                  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),plot.background=element_blank()) +
            geom_density_2d(data=projectionDF[projectionDF$Type == "organoid",], mapping=aes(x=X1,y=X2),color="#782867",n=300,h=3)) + NoLegend() + NoAxes()
    dev.off()
    
  } else if (!is.integer0(grep("MID", id))){
    
    ## MID
    pdf(paste0(plotDir,id,"-Projection-UMAP.pdf"), width = 12, height = 12, useDingbats = FALSE)
    print(ggplot(projectionDF, aes(X1,X2,color=Type)) + 
            geom_point() +
            theme_bw() +
            xlab("UMAP Dimension 1") + 
            ylab("UMAP Dimension 2") +
            scale_color_manual(values=c("reference"="lightgrey", "organoid" = "darkblue")) +
            theme_classic() +
            theme(axis.line=element_blank(),axis.text.x=element_blank(),
                  axis.text.y=element_blank(),axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),legend.position="none",
                  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),plot.background=element_blank()) +
            geom_density_2d(data=projectionDF[projectionDF$Type == "organoid",], mapping=aes(x=X1,y=X2),color="darkblue",n=300,h=3)) + NoLegend() + NoAxes()
    dev.off()
    
  } else{
    
    ## WT
    pdf(paste0(plotDir,id,"-Projection-UMAP.pdf"), width = 12, height = 12, useDingbats = FALSE)
    print(ggplot(projectionDF, aes(X1,X2,color=Type)) + 
            geom_point() +
            theme_bw() +
            xlab("UMAP Dimension 1") + 
            ylab("UMAP Dimension 2") +
            scale_color_manual(values=c("reference"="lightgrey", "organoid" = "#4F9E4C")) +
            theme_classic() +
            theme(axis.line=element_blank(),axis.text.x=element_blank(),
                  axis.text.y=element_blank(),axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),legend.position="none",
                  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),plot.background=element_blank()) +
            geom_density_2d(data=projectionDF[projectionDF$Type == "organoid",], mapping=aes(x=X1,y=X2),color="#4F9E4C",n=300,h=3)) + NoLegend() + NoAxes()
    dev.off()
  }
  
  #############################################################
  #PERFORM DIFFERENTIAL GENE EXPRESSION ANALYSIS IN LSI SPACE
  #############################################################
  
  #LSI-SVD
  svdReference <- as.data.frame(lsiReference$matSVD)
  svdDisease <- as.data.frame(as.matrix(lsiProjection))
  
  #Differential Seed
  set.seed(1)
  
  #Cells that we are testing of disease
  idxDisease <- rownames(projectionDF)[projectionDF$CellType=="organoid"]
  
  #If the number of cells exceeds the max downsample to max
  if(length(idxDisease) > nMax){
    idxDisease <- sample(idxDisease, nMax)
  }
  
  #If the number of cells is greater than 1 continue
  stopifnot(length(idxDisease) > 1)
  
  #KNN Nearest Neighbor using FNN
  knnDisease <- FNN::get.knnx(
    data = svdReference,
    query = svdDisease[idxDisease, ], #Subset by idxDisease 
    k = input_knn)
  
  # nn.index - an n x k matrix for the nearest neighbor indice.
  # make a nearest neighbor dataframe
  nn.df <- data.frame(nn_index = as.vector(knnDisease$nn.index), Cell_barcodes = str_split_fixed(rownames(svdReference[as.vector(knnDisease$nn.index),]), "-", 2)[,1])
  nn.df <- nn.df[order(nn.df$nn_index),]
  nn.df$Cell_barcodes <- paste0(nn.df$Cell_barcode, "-1")
  
  # merge the nn.df with the projectionDF from before to get information about cell identity etc.
  full.nn.df <- merge(nn.df, projectionDF, by = "Cell_barcodes", all.x = T)
  
  # nn overview of cells for the sample
  nn.frequency <- table(full.nn.df$CellType) %>% bind_rows
  nn.table <- rBind.fill(nn.table, as.data.frame(nn.frequency))
  rownames(nn.table)[i+1] <- id
  
  #Determine the minimum KNN where reference cells are less than 1.25x disease cells
  i <- 0
  uniqueIdx <- unique(as.vector(knnDisease$nn.index))
  while(length(uniqueIdx) > 1.25 * length(idxDisease)){
    i <- i + 1
    uniqueIdx <- unique(as.vector(knnDisease$nn.index[,seq_len(input_knn-i)]))
  }
  
  #Reference cells for testing
  idxReference <- rownames(svdReference)[uniqueIdx]
  
  #If there are more healthy cells downsample healthy cells
  #If there are more disease cells downasmple disease cells
  if(length(idxReference) > length(idxDisease)){
    idxReference <- sample(idxReference, length(idxDisease))
  }else{
    idxDisease <- sample(idxDisease, length(idxReference))
  }
  message(sprintf("nDisease = %s\nnHealthy = %s", length(idxDisease), length(idxReference)))
  
  #Disease and Reference Matrix
  matHealthy <- assay(se[,idxReference])
  matDisease <- assay(seDisease[,idxDisease])
  
  #Normalize to scaleTo
  matNormDisease <- t(t(matDisease)/Matrix::colSums(matDisease)) * scaleTo
  matNormHealthy <- t(t(matHealthy)/Matrix::colSums(matHealthy)) * scaleTo
  
  #T-Test Comparisons
  dfTT <- sparseMatTTest(matNormDisease, matNormHealthy)
  dfTT$feature <- rownames(matNormDisease)
  dfTT$log2Mean <- log2(rowMeans(cbind(dfTT$mean1, dfTT$mean2)) + 10^-4)
  dfTT$log2FC <- log2((dfTT$mean1 + 10^-4)/(dfTT$mean2 + 10^-4))
  
  plotDiff <- data.frame(row.names=row.names(dfTT),log2Mean=dfTT$log2Mean,log2FC=dfTT$log2FC,FDR=dfTT$fdr)
  plotDiff <- plotDiff[complete.cases(plotDiff),]
  plotDiff$type <- "not-differential"
  plotDiff$type[plotDiff$log2FC > 0.5 & plotDiff$FDR < 0.01] <- "up-regulated"
  plotDiff$type[plotDiff$log2FC < -0.5 & plotDiff$FDR < 0.01] <- "do-regulated"
  
  
  pdf(paste0(plotDir,id,"-Differential-MA-Plot.pdf"), width = 8, height = 6, useDingbats = FALSE)
  print(ggplot(plotDiff, aes(log2Mean,log2FC,color=type)) + 
          geom_point(size=0.5) +
          theme_bw() +
          xlab("log2 Mean") + 
          ylab("log2 Fold Change") +
          scale_color_manual(values=c("not-differential"="lightgrey", "do-regulated"="dodgerblue3", "up-regulated"="firebrick3")))
  dev.off()
  
  #Save Output
  readr::write_tsv(dfTT, paste0(plotDir,id,"-Differential-Results.tsv"))
  
}

# write the KNN table to file
write.table(nn.table, paste0(sample.tmp, "/", sample.tmp, "_KNN_celltype_table.txt"), row.names = T, col.names = T, quote = F, sep = "\t")

################################################################
#       CREATE A PLOT FOR THE DRIFT QUANTIFICATION
################################################################
# nn.table <- read.delim(paste0(sample.tmp, "/", sample.tmp, "_KNN_celltype_table.txt"), header = T)

colnames(nn.table) <- c("Chief Cells", "PMCs", "Enterocytes", "Parietal Cells", "Enteroendocrine Cells", "Goblet Cells", "MSCs", "GMCs", "Mucosal-like malignant Cells", 
                        "Non-Mucosal-like malignant Cells", "Proliferative Cells", "Paneth Cells")

# calculate percentages of nearest neighbors 
# then compare for each sample - early late
# is the frequency of cells mapping to cancer increasing
nn.table <- nn.table[-1,]

# calculate frequencies instead of number of nearest neighbors
sample.cells <- rowSums(nn.table)
freq.table <- round(nn.table/sample.cells, 3)*100

# change D2 WT
freq.table <- rbind(freq.table, freq.table[grep("D2_WT", rownames(freq.table)),], freq.table[grep("D2_WT", rownames(freq.table)),])
rownames(freq.table)[nrow(freq.table)-1] <- "D2_C1_WT"
rownames(freq.table)[nrow(freq.table)] <- "D2_C2_WT"
rownames(freq.table)[grep("D2_WT", rownames(freq.table))] <- "D2_C3_WT"

# change D1 WT
freq.table <- rbind(freq.table, freq.table[grep("D1_WT", rownames(freq.table)),], freq.table[grep("D1_WT", rownames(freq.table)),])
rownames(freq.table)[nrow(freq.table)-1]  <- "D1_C1_WT"
rownames(freq.table)[nrow(freq.table)]  <- "D1_C2_WT"
rownames(freq.table)[grep("D1_WT", rownames(freq.table))] <- "D1_C3_WT"

# change D3 WT
freq.table <- rbind(freq.table, freq.table[grep("D3_WT", rownames(freq.table)),], freq.table[grep("D3_WT", rownames(freq.table)),])
rownames(freq.table)[nrow(freq.table)-1] <- "D3_C1_WT"
rownames(freq.table)[nrow(freq.table)] <- "D3_C2_WT"
rownames(freq.table)[grep("D3_WT", rownames(freq.table))] <- "4230_C3_WT"

# add sample information
freq.table$sample_id <- as.character(rownames(freq.table))

# melt the freq table to plot it 
melt.freq.table <- melt(as.data.frame(freq.table), id.var = "sample_id")

# add timepoint information
melt.freq.table$timepoint <- "WT"
melt.freq.table[grep("EARLY", melt.freq.table$sample_id),"timepoint"] <- "EARLY"
melt.freq.table[grep("MID", melt.freq.table$sample_id),"timepoint"] <- "MID"
melt.freq.table[grep("LATE", melt.freq.table$sample_id),"timepoint"] <- "LATE"
melt.freq.table$timepoint <- factor(melt.freq.table$timepoint, levels = c("WT", "EARLY", "MID", "LATE"))

# add clone information
melt.freq.table$clone_id <- paste0(str_split_fixed(melt.freq.table$sample_id, "_", 3)[,1], "_", str_split_fixed(melt.freq.table$sample_id, "_", 3)[,2])

# replace the cell names that are not correctly shown
melt.freq.table$variable <- as.character(melt.freq.table$variable)

# rename some cell types
melt.freq.table[grep("PMCs", melt.freq.table$variable),"variable"] <- "Pit Mucosal Cells"
melt.freq.table[grep("GMCs", melt.freq.table$variable),"variable"] <- "Gland Mucosal Cells"
melt.freq.table[grep("MSCs", melt.freq.table$variable),"variable"] <- "Mucosal Stem Cells"

# order the cell types
melt.freq.table$variable <- factor(melt.freq.table$variable, levels = c("Mucosal-like malignant Cells", "Non-Mucosal-like malignant Cells",
                                                                        "Pit Mucosal Cells", "Gland Mucosal Cells", 
                                                                        "Chief Cells", "Mucosal Stem Cells",
                                                                        "Enteroendocrine Cells", "Endocrine Cells",
                                                                        "Enterocytes", "Goblet Cells", 
                                                                        "Paneth Cells", "Parietal Cells", "Proliferative Cells"))

# check maximum 
max.celltype.contribution <- aggregate(melt.freq.table$value, by = list(melt.freq.table$variable), max)
valid.celltypes <- as.character(max.celltype.contribution[max.celltype.contribution$x > 15, "Group.1"])

# remove cell types with a maximum contribution of 5%
melt.freq.table <- melt.freq.table[melt.freq.table$variable %in% valid.celltypes,]
melt.freq.table$donor <- str_split_fixed(melt.freq.table$clone_id, "C", 2)[,1]
melt.freq.table$clone <- substring(melt.freq.table$clone_id, 3, 4)

color = c("D1"="#00B050","D2"="#843C0C","D3"="#7C7C7C") ### COLORFUL

# plot data
ggplot(data = melt.freq.table, aes(x = timepoint, y = value, group = clone_id, shape = clone)) +
  geom_line(aes(color = donor), size = 1) +
  geom_point(aes(color = donor), size = 3) +
  scale_color_manual(values=color) +
  facet_wrap(~ variable, ncol = 4) +
  theme_bw() + 
  guides(alpha = FALSE) +
  labs(title = "K-Nearest Neighbor Cell Type Quantification", x = "Time", y = "Cell type frequency [%]", color = "Clone IDs") +
  theme(strip.text = element_text(face="bold", size=18, colour = "black",),
        strip.background = element_rect(fill="lightgrey", colour="black", size=1), 
        axis.text = element_text(colour = "black", size = 14, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 14, face = "bold" ),
        axis.title = element_text(colour = "black", size = 20, face = "bold" ),
        plot.title = element_text(colour = "black", size = 20, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 14, face = "bold",),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "bottom") +
  geom_text_repel(data = melt.freq.table %>% filter(timepoint == "LATE") %>% group_by(clone_id) %>% filter(value == max(value)),
                  aes(label = paste0(clone_id, " - ", value, "%")) ,
                  hjust = -.35,
                  nudge_x = .5,
                  direction = "y",
                  fontface = "bold",
                  size = 5) +
  guides(color=guide_legend( override.aes = list(size = 3)))
ggsave(paste0(sample.tmp, "/AllPatients_LSI_quantification_corrected_", sample.tmp ,".pdf"), width = 18, height = 9)
