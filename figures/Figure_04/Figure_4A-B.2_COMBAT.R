#############################################################################################################################
##                                                                                                                      
##  Remove organoid batch effects using the sva package (Johnson et al., 2007)                                                                                          
##                                                                                                                      
##  Date: 20 April 2020                                                                                                                    
##  
##  Author: Moritz Przybilla                                                                                                                    
##                                                                                                                      
############################################################################################################################

# clear workspace
rm(list=ls())
set.seed(14) # set the seed of random number generator 

# package dependencies, which have to be installed are checked and installed if nnot available
list.of.packages <- c("tidyverse", "patchwork", "Seurat", "Matrix", "sva", "bladderbatch", "pamr", "limma")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

#####################################################################################
# FUNCTIONS
#####################################################################################

# merge the matrices together
check_and_merge <- function(sparse_matrix_list){
  if(length(sparse_matrix_list) < 2){
    return(sparse_matrix_list[[1]])
  }
  combtab <- combn(x = seq(length(sparse_matrix_list)),m = 2)
  
  genes_in_common <- list()
  for(i in seq(ncol(combtab))){
    genes_in_common[[i]] <- intersect(rownames(sparse_matrix_list[[combtab[1,i]]]), rownames(sparse_matrix_list[[combtab[2,i]]]))
  }
  
  cg <- genes_in_common[[1]]
  if(length(genes_in_common) > 1 ){
    for(i in 2:length(genes_in_common)){
      cg <- intersect(cg,genes_in_common[[i]])
    }
  }
  
  filter_genes <- function(x, genes){
    x <- x[genes,]
    return(x)
  }
  
  sparse_matrix_list <- lapply(sparse_matrix_list, filter_genes, cg)
  
  check.cols <- c()
  for(i in seq(ncol(combtab))){
    check.cols <- c(check.cols, identical(rownames(sparse_matrix_list[[combtab[1,i]]]), rownames(sparse_matrix_list[[combtab[2,i]]])) )
  }
  if(all(check.cols)){
    outmat <- do.call(cbind, sparse_matrix_list)
    return(outmat)
  } else {
    return(NULL)
  }
}

#####################################################################################
# READ IN DATA OF INTEREST
#####################################################################################

# read in the epithelial Zhang dataset
Sathe.epithelial.seurat.obj <- readRDS("/labs/ccurtis2/mjprzy/scRNA_analysis/hashECB_data_freeze/merged_data_SatheEtal/merged_data_SatheEtal_Epithelial_CellType_seurat_obj.rds")
Sathe.matrix <- Sathe.epithelial.seurat.obj@assays$RNA@counts

# show dimensions
print(dim(Sathe.matrix))

## GET THE INDIVIDUAL PATIENT ORGANOID DATA AND COMBINE THE SEURAT OBJECTS
EML.seurat.obj.big <- readRDS("/labs/ccurtis2/mjprzy/scRNA_analysis/hashECB_data_freeze/EML_Seq_freeze/EML_Seq_freeze_seurat_obj.rds")
EML.matrix <- EML.seurat.obj.big@assays$RNA@counts

# show dimensions
print(dim(EML.matrix))

# make a metadata sheet including cell_id outcome batch and wether its organoid or tissue with cell barcodes as rownames
EML.organoid.metadata <- data.frame(sample_id = EML.seurat.obj.big@meta.data[,2])
rownames(EML.organoid.metadata) <- rownames(EML.seurat.obj.big@meta.data)
EML.organoid.metadata$cell <- c(1:nrow(EML.organoid.metadata))
EML.organoid.metadata$batch <- 1
EML.organoid.metadata$type <- "Organoid"

# for the Sathe dataset
Sathe.epithelial.metadata <- data.frame(sample_id = Sathe.epithelial.seurat.obj@meta.data[,1])
rownames(Sathe.epithelial.metadata) <- rownames(Sathe.epithelial.seurat.obj@meta.data)
Sathe.epithelial.metadata$cell <- c((nrow(EML.organoid.metadata)+1):(nrow(EML.organoid.metadata)+nrow(Sathe.epithelial.metadata)))
Sathe.epithelial.metadata$batch <- 0
Sathe.epithelial.metadata$type <- "Tissue"

# bind info together
s.info <- rbind(EML.organoid.metadata, Sathe.epithelial.metadata)

# define output directory
o.dir <- "/labs/ccurtis2/mjprzy/scRNA_analysis/combat/freeze"
dir.create(o.dir)

# define sample name
sample.tmp <- "EML_Seq_freeze"

# create o.dir for sample
dir.create(paste0(o.dir, "/", sample.tmp))
setwd(paste0(o.dir, "/", sample.tmp))

# make a plot directory
dir.create("plots")

###########################################################################
#                       Getting started with COMBAT
###########################################################################
#' The COMBAT function adjusts for known batches byusing an empirical Bayesian framework.
#' In theory, COMBAT should remove signals that distinguish organoids from cancer tissues
#' as a group, which is essentially the signal coming from the tumor micro-
#' environment. This might nnormally dominate the clustering, but should be removed after
#' use of COMBAT. 
#' 
#' The sva package assumes the availability of a matrix with features x samples. 
#' Furthermore, it assumes that there are two types of varaibles that are being
#' considered: (1) adjustment variables (i.e. age of the patients) and 
#' (2) variables of interest (i.e. indicator of cancer vs. control). 
#' 
#' Two model matrices must be made:
#' 1. "Null model"-matrix:
#'    Model matrix that includes terms for all of the adjustment variables but not the
#'    variables of interest.
#' 2. "Full model"-matrix:
#'    Model matrix taht includes terms for both the adjustment variables and the 
#'    variables of interest.

print('Batch correcting with ComBat')

# use a known batch variable in your dataset
# the batch variable is a vector identical to the number of columns in your matrix
batch.use <- s.info$batch # batch.use = vector of batch labels for each cell

# merge the matrices from Organoids and Zhang together
matrix.list <- list(EML.matrix, Sathe.matrix)
combined.matrix <- check_and_merge(matrix.list)

###########################################################################
#           Create a model matrix for the adjustment varibles
###########################################################################
# use the model.matrix function to set these up
modcombat <- model.matrix(~1, data = s.info)

#' Adjustment variables will be treated as given to the COMBAT function.

###########################################################################
#                       Run the COMBAT function
###########################################################################
# Apply the COMBAT function to the data, using parametric empirical Bayesian adjustments
# bc - batch correction
combat.bc.data <- ComBat(dat=as.matrix(combined.matrix), batch.use, par.prior=TRUE, prior.plots=FALSE)

#' This returns an expression matrix, with the same dimensions as your original dataset.

# # write corrected data to file
# write.table(combat.bc.data, file=paste0(sample.tmp, '_bc_data.txt'), sep='\t', quote = F, col.names = T, row.names = T)

# now, get the organoid cells only
EL.combat.bc.data <- combat.bc.data[, colnames(combat.bc.data) %in% colnames(EML.matrix)]

# write corrected data to file
write.table(EL.combat.bc.data, file=paste0(sample.tmp, '_organoid_bc_data.txt'), sep='\t', quote = F, col.names = T, row.names = T)


