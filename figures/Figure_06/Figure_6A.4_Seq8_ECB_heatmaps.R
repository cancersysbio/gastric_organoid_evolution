################################################################################################################################################
##                                                                                                                      
##  Plot subset of gene expression matrix based on the source code from honeybadger (credit to Jean Fan) for Sequencing 8                                                               
##                                                                                                                      
##  Date: 11 April 2020                                                                                                                    
##  
##  Author: Moritz Przybilla                                                                                                    
##                                                                                                                      
##                                                                                                                      
################################################################################################################################################
# clear workspace beforehand
rm(list = ls())
set.seed(2020)

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("biomaRt", "Matrix", "stringr", "GenomicRanges", "scater")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

#####################################################################################
#                                  load data
#####################################################################################
# read in normalized matrix
norm.gExp.matrix <- as.matrix(read.table("/labs/ccurtis2/mjprzy/infercnv_gastric/freeze/Sequencing8_C5T2R2_cr3_WT/results/infercnv.median_filtered.observations.txt", header = T, sep = " "))
colnames(norm.gExp.matrix) <- str_split_fixed(colnames(norm.gExp.matrix), "\\.", 2)[,1]

# read in annotation file
annotation.file <- read.table("/labs/ccurtis2/mjprzy/infercnv_gastric/freeze/Sequencing8_C5T2R2_cr3_WT/cellAnnotations_Sequencing8_C5T2R2_cr3_WT.txt", header = F, sep = "\t", stringsAsFactors = F)
annotation.file$V1 <- str_split_fixed(annotation.file$V1, "-", 2)[,1]
annotation.file <- annotation.file[annotation.file$V2 != "normal",]

# read in gene_ordering file
gene.ordering.file <- read.table("/labs/ccurtis2/mjprzy/infercnv_gastric/freeze/Sequencing8_C5T2R2_cr3_WT/gene_ordering_file.txt", header = F, col.names = c("hgnc_symbol", "chr", "start", "end"), sep = "\t")
gene.ordering.file$coordinates <- paste(paste0("chr", gene.ordering.file$chr), gene.ordering.file$start, gene.ordering.file$end, sep = ":")
gene.ordering.file$chr <- paste0("chr", gene.ordering.file$chr)
rownames(gene.ordering.file) <- gene.ordering.file$hgnc_symbol

# read in metadata
seurat.obj <- readRDS("/labs/ccurtis2/mjprzy/scRNA_analysis/hashECB_data_freeze/Sequencing8_C5T2R2_cr3/Sequencing8_C5T2R2_cr3_seurat_obj.rds")
metadata <- seurat.obj@meta.data
metadata <- metadata[,c(1,5,6,7)]

# ECB clones
ECB.clones <- table(metadata$ECB_RG) > 10
valid.clones <- names(ECB.clones[ECB.clones == T])

#####################################################################################
#                   Generate input for the image function to plot
#####################################################################################
# center the infercnv output to 0
gexp.norm <- norm.gExp.matrix -1

# only take genes which are present in the matrix
genes <- makeGRangesFromDataFrame(gene.ordering.file)
genes <- genes[rownames(gexp.norm)]

# generate new object genes.interest with genes and respective coordinates
genes.interest <- as.data.frame(genes)
rownames(genes.interest) <- names(genes)
mat <- gexp.norm

# check for each chromosome if there are genes present
gExp.array <- tapply(1:nrow(genes.interest),as.factor(genes.interest$seqnames),function(ii) {
  na.omit(mat[rownames(genes.interest)[ii[order((genes.interest[ii,]$start+genes.interest[ii,]$end)/2,decreasing=F)]],,drop=FALSE])
})

# subset to chromosome 3,4,9,13 and 20 (chr of interest in general)
chrs <- paste0("chr", 1:22)
gExp.array <- gExp.array[chrs]

# we want to scale the output to the respective amount of genes and their widths
widths <- sapply(gExp.array, nrow); widths <- widths/max(widths)*100

# generate the layout accordingly
plot.layout <- layout(matrix(seq(1,length(gExp.array)),1,length(gExp.array),byrow=TRUE), widths=widths)

# calculate average distance
groups <- valid.clones
chr.list <- list()
ordered_names <- c()

for (j in 1:length(groups)){
  
  # which group?
  print(groups[j])
  
  # get the barcodes from the respective group
  barcodes <- str_split_fixed(metadata[grep(paste("^",groups[j],"$", sep=""), metadata$ECB_RG),"Cell_Barcode"], "-1", 2)[,1]
  
  # iterate over all chrs and subset the matrix
  for (i in 1:length(gExp.array)){
    
    # matrix 
    matrix <- gExp.array[[i]]
    matrix <- matrix[,colnames(matrix) %in% barcodes]
    
    chr.list[[paste0("chr",i)]] <- matrix
  }
  
  # calculate the col averages here for the subset
  avgd <- do.call(rbind, lapply(names(chr.list),function(nam) {
    d <- chr.list[[nam]]
    d <- colMeans(d)
    d
  }))
  
  # Order cells with hierarchical clustering
  dist.centered.matrix <- dist(t(avgd), method = "euclidean")
  hc <- hclust(dist.centered.matrix, method = "ward.D2")
  
  # make a vector with the right order of barcodes per group
  ordered_names <- c(ordered_names, hc$labels[hc$order])
}

# set parameters for the image
adapted.gExp.list <- gExp.array
pcol <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))(256) # color panel
zlim <- c(-0.2, 0.2) # upper and lower limit

# limit the gExp to maximums according to zlim
limit.gExp.list <- lapply(names(adapted.gExp.list),function(nam) {
  d <- adapted.gExp.list[[nam]]
  d[d< zlim[1]] <- zlim[1]; d[d>zlim[2]] <- zlim[2];
  return(d)
})

#####################################################################################
#                   Generate ECB colour coding in the right order
#####################################################################################

# first, remove the malignant to free the barcode id
annotation.file$ECB_RG <- as.numeric(str_split_fixed(annotation.file$V2, "malignant_", 2)[,2])

# make an dataframe with an index for each ECB_RG
unique.ECB_RG <- data.frame(c(unique(annotation.file$ECB_RG)))
unique.ECB_RG$index <- rownames(unique.ECB_RG)
colnames(unique.ECB_RG) <- c("ECB_RG","index")

# the index is important for the order so we merge it with the seurat metadata
ECB.hc.metadata <- merge(metadata, unique.ECB_RG, by="ECB_RG")

# and make the ECBs as factors to order the dataframe according to them, as well as the hc order
ECB.order <- c(0:10000)
ECB.hc.metadata$ECB_RG <- factor(ECB.hc.metadata$ECB_RG, levels = ECB.order)

# remove -1 to make it consistent
ECB.hc.metadata$Cell_Barcode <- str_split_fixed(ECB.hc.metadata$Cell_Barcode, "-1", 2)[,1] 

# and order it finally to the ECB and within each ECB hierarchically
final.metadata <- ECB.hc.metadata

# make a df with the index to color it accordingly
ECB.col.bar <- data.frame(as.numeric(final.metadata$index),as.numeric(final.metadata$index))
colnames(ECB.col.bar)<-c("rg1","rg2")

# get the coloring according to RGs
cols <- c(as.character(unique(final.metadata$RG_color)))

#####################################################################################
#                                 Make the plot
#####################################################################################

png("/labs/ccurtis2/mjprzy/infercnv_gastric/freeze/Sequencing8_C5T2R2_cr3_WT/Sequencing8_C5T2R2_cr3_WT_w_ECB.png" , width = 26, height = 12, units = 'in', res = 600)
par(mfrow=c(1,23), mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.8)
image(t(ECB.col.bar),xlab="",ylab="", axes=FALSE, col=cols)
## plot chromosomes
box()
for (i in 1:length(limit.gExp.list)){
  message(paste0("chr", i))
  d <- limit.gExp.list[[i]]
  d <- d[, ordered_names] 
  # image the respective chr 
  image(seq_len(nrow(d)), seq_len(ncol(d)), d, col=pcol, zlim=zlim, xlab="", ylab="", axes=F, main=paste0("chr", i), useRaster = T)
  box()
}
dev.off()

png("/labs/ccurtis2/mjprzy/infercnv_gastric/freeze/Sequencing8_C5T2R2_cr3_WT/Sequencing8_C5T2R2_cr3_WT_w_ECB_legend.png" , width = 13, height = 4, units = 'in', res = 600, type = "cairo")
par(mfrow=c(1,24), mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.5)
image.plot(legend.only=TRUE, zlim=zlim, col = pcol, legend.shrink = 0.4, legend.width = 10, smallplot=  c(.65, .9, .5, .8)) # first argument of small plot = width of legend, third is making the length, the higher the smaller
dev.off()

