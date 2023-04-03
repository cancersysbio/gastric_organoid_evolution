#!/usr/bin/env Rscript

## Author: Aziz Khan - <azizk@stanford.edu>

## Version 2.0
# Quantitative DNA sequencing for chromosomal aberrations. The genome is divided into non-overlapping fixed-sized bins, 
# number of sequence reads in each counted, adjusted with a simultaneous two-dimensional loess correction for sequence mappability 
# and GC content, and filtered to remove spurious regions in the genome.
# Downstream steps of segmentation and calling are also implemented via packages DNAcopy and CGHcall, respectively.

suppressPackageStartupMessages({
library("optparse")
library(tidyverse)
library(Biobase) 
library("QDNAseq")
library(CGHbase)
future::plan("multiprocess")
library(data.table)
library(ACE)
library(png)
library(grid)
library(ggplot2)
library(gridExtra)
})


option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Input BAM file or path", metavar="character"),
	make_option(c("-o", "--out"), type="character", default="./", 
              help="output path name [default= %default]", metavar="character"),
  make_option(c("-n", "--name"), type="character", default="sample1", 
              help="sample bam file name/id [default= %default]", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default="hg38", 
              help="Genome (hg19, hg38) [default= %default]", metavar="character"),
  make_option(c("-m", "--minmapq"), type="integer", default=37, 
              help="Minimum quality score [default= %default]", metavar="integer"),
  make_option(c("-r", "--refit"), type="logical", default=TRUE, 
              help="Refit the calls using ACE package [default= %default]", metavar="logical"),
   make_option(c("-d", "--duplicated"), type="logical", default=FALSE, 
              help="Are the reads dublicated? [default= %default]", metavar="logical"),
	make_option(c("-b", "--bin"), type="integer", default=50,
              help="bin size (kb). Use 5, 10, 15, 30, 50, 100, 200, 500, 1000 [default= %default]", metavar="integer")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input bam file).n", call.=FALSE)
}

##set working dir
setwd(opt$out)
opt$out <- paste0(opt$out,"/")
sample_prefex <- paste0(opt$name,"_",opt$bin)
outdir_prefex <- paste0(opt$out,sample_prefex)

sink(paste0(outdir_prefex,"_log.txt"))

## load the bin annotations
if (opt$genome == 'hg38'){
  library("QDNAseq.hg38")
  bins <- getBinAnnotations(binSize=opt$bin, genome="hg38")
 
  }else{
    library("QDNAseq.hg19")
    bins <- getBinAnnotations(binSize=opt$bin, genome="hg19")
}

probability_matrix <- function(cgh) {
    probs = matrix(0, dim(cgh)[1], dim(cgh)[2])
    colnames(probs) = colnames(cgh)
    rownames(probs) = rownames(cgh)

    for (i in 1:dim(cgh)[1]){
        for (j in 1:dim(cgh)[2]){
            if (calls(cgh)[i,j] == 2){
                probs[i,j] = probamp(cgh)[i,j]
            } else if (calls(cgh)[i,j] == 1){
                probs[i,j] = probgain(cgh)[i,j]
            } else if (calls(cgh)[i,j] == -1){
                probs[i,j] = probloss(cgh)[i,j]
            } else if (calls(cgh)[i,j] == -2){
                probs[i,j] = probdloss(cgh)[i,j]
            }
        }
    }

    probs
}

getSegTable <- function(x)
{
  dat<-x
  sn<-assayDataElement(dat,"segmented")
  fd <- fData(dat)
  fd$use -> use
  fdfiltfull<-fd[use,]
  sn<-sn[use,]
  segTable<-c()
  for(c in unique(fdfiltfull$chromosome))
  {
    snfilt<-sn[fdfiltfull$chromosome==c]
    fdfilt<-fdfiltfull[fdfiltfull$chromosome==c,]
    sn.rle<-rle(snfilt)
    starts <- cumsum(c(1, sn.rle$lengths[-length(sn.rle$lengths)]))
    ends <- cumsum(sn.rle$lengths)
    lapply(1:length(sn.rle$lengths), function(s) {
      from <- fdfilt$start[starts[s]]
      to <- fdfilt$end[ends[s]]
      segValue <- sn.rle$value[s]
      c(fdfilt$chromosome[starts[s]], from, to, segValue)
    }) -> segtmp
    segTableRaw <- data.frame(matrix(unlist(segtmp), ncol=4, byrow=T),stringsAsFactors=F)
    segTable<-rbind(segTable,segTableRaw)
  }
  colnames(segTable) <- c("chrom", "start", "end", "segVal")
  return(segTable)
}

if (file_test("-f", opt$file)){
    bams <- c(opt$file)
  }else{
    bams <- paste0(opt$file, list.files(opt$file, pattern=".bam", recursive=F))
}


#if (!file.exists(paste0(outdir_prefex,"_copyNumbersSegmented.rds"))){

readCounts <- binReadCounts(bins,  bamfiles=bams, bamnames=c(opt$name), minMapq=opt$minmapq, isDuplicate=opt$duplicated)
#saveRDS(readCounts, file=paste0(outdir_prefex,"_readCounts.Rdata"))

pdf(paste0(outdir_prefex,"_rawprofile.pdf"), width=20, height=8, useDingbats=F)
## Plot the readcounts with filtered reads highlighted.
plot(readCounts, logTransform=TRUE)
highlightFilters(readCounts, logTransform=FALSE,residual=TRUE, blacklist=TRUE)
dev.off()

##Apply QDNAseq filters.
readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)
# fix for issue with copy number segmentation arising from zero count bins
# plot median read counts as a function of GC content and mappability as an isobar plot
png(paste(outdir_prefex,"_isobar.png", sep=''))
isobarPlot(readCountsFiltered)
dev.off()

##Calculate CG correction.
readCountsFiltered <- estimateCorrection(readCountsFiltered)
png(paste(outdir_prefex,"_noise.png",sep=''))
noisePlot(readCountsFiltered)
# noise plot showing relationship between the observed standard deviation in the data
# and its read depth
dev.off()

##Apply GC correction.
copyNumbers <- correctBins(readCountsFiltered)
#saveRDS(copyNumbers, file=paste0(outdir_prefex,"_readCountsFilteredCorrected.Rdata"))

##Normalise and smooth outliers.
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

pdf(paste0(outdir_prefex,"_copyNumbersSmooth.pdf"), width=20, height=8, useDingbats=F)
##Plot the smoothed copy-number.
plot(copyNumbersSmooth)
dev.off()


# output raw and fitted read counts
features <- fData(copyNumbersSmooth) %>%
    as.data.frame %>%
    rownames_to_column(var = "location") %>%
    transmute(location, chrom = chromosome, start = as.integer(start), end = as.integer(end))
  
rawReadCounts <- assayData(copyNumbersSmooth)$counts %>%
    as.data.frame %>%
    rownames_to_column(var = "location")
  features %>%
    dplyr::left_join(rawReadCounts, by = "location") %>%
    write_tsv(paste(outdir_prefex,"_rawReadCounts.txt",sep=''))
  
fittedReadCounts <- assayData(copyNumbersSmooth)$fit %>%
    as.data.frame %>% rownames_to_column(var = "location") %>%
    mutate_if(is.numeric, list(~round(., digits = 3)))
  features %>%
    dplyr::left_join(fittedReadCounts, by = "location") %>%
    write_tsv(paste(outdir_prefex,"_fittedReadCounts.txt",sep=''))

## Segment the copy-number profile. 
copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")

copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

##Plot the segmented profile.

pdf(paste0(outdir_prefex,"_CN_segmented.pdf"), width=15, height=5, useDingbats=F)
##PDF plot
plot(copyNumbersSegmented)
dev.off()


##png plot
png(paste0(outdir_prefex,"_CN_segmented.png"), width = 1200, height = 300)
plot(copyNumbersSegmented)
dev.off()

pData(copyNumbersSegmented) %>%
    rownames_to_column(var = "sample") %>%
    dplyr::select(id = name, everything()) %>%
    write_tsv(paste0(outdir_prefex,"_SegmentedSummary.txt"))

saveRDS(copyNumbersSegmented, file=paste0(outdir_prefex,"_copyNumbersSegmented.rds"))
#}else
#{
#  copyNumbersSegmented <- readRDS(paste0(outdir_prefex,"_copyNumbersSegmented.rds"))
#}

## Run ACE

ploidies = c(2,3,4)
ACE_output_dir <- paste0(outdir_prefex,"_ACE_fits/")
#if(opt$refit == TRUE){
if (!dir.exists(ACE_output_dir)){
 system(paste0("rm -r ",outdir_prefex,"_ACE_fits"))
 runACE(opt$out, outputdir=ACE_output_dir, filetype = 'rds', genome = opt$genome, 
  binsizes=opt$bin, ploidies = ploidies, imagetype = 'png', method = 'RMSE', penalty = 0, 
  cap = 12, bottom = 0, trncname = ".recal", printsummaries = TRUE, autopick = TRUE)
}
#}

ace_fits <- data.frame()
error_plots <- c()
for (ploidy in ploidies){
  ploidy <- paste0(ploidy,"N")
  ace_fits <-  rbind(ace_fits, read.delim(paste0(ACE_output_dir,sample_prefex,"_copyNumbersSegmented/",ploidy,"/fitpicker_",ploidy,".tsv")))
  error_plots <- append(error_plots, paste0(ACE_output_dir,sample_prefex,"_copyNumbersSegmented/",ploidy,"/summary_errors.png"))
}
write.table(ace_fits, paste0(ACE_output_dir,sample_prefex,"_ACE_likely_fits.txt"), sep='\t', quote=FALSE, row.names=FALSE)

## merge the ACE likelyfit error plots into one
plots <- lapply(ll <- error_plots ,function(x){
  img <- as.raster(readPNG(x))
  rasterGrob(img, interpolate = FALSE)
})
ggsave(paste0(ACE_output_dir,sample_prefex,"_ACE_summary_errors.png"), width=12, height=12, marrangeGrob(grobs = plots, nrow=3, ncol=1, top=NULL))


cellularity <- max(ace_fits$likely_fit)


##QDNAseq use CGHcall for CN calls and below are default arguments
#CGHcall(inputSegmented, prior = "auto", nclass = 5, organism = "human", cellularity=1, robustsig="yes", nsegfit=3000, maxnumseg=100, minlsforfit=0.5, build="GRCh37", ncpus=1)

## Call aberrations from segmented copy number data

if (opt$genome == 'hg38'){
  copyNumbersCalled <- callBins(copyNumbersSegmented, cellularity=cellularity, build="GRCh38")
  }else{
    copyNumbersCalled <- callBins(copyNumbersSegmented, cellularity=cellularity, build="GRCh37")
}

exportBins(copyNumbersCalled, file=paste0(outdir_prefex,"_CNA_called.vcf"), format="vcf")
exportBins(copyNumbersCalled, file=paste0(outdir_prefex,"_CNA_called.seg"), format="seg")


## Plot final profile.
pdf(paste0(outdir_prefex,"_CN_called.pdf"), width=15, height=5, useDingbats=F)
plot(copyNumbersCalled)
dev.off()

pdf(paste0(outdir_prefex,"_frequency_plot.pdf"), width=15, height=5, useDingbats=F)
frequencyPlot(copyNumbersCalled)
dev.off()

cgh <- makeCgh(copyNumbersCalled)
write.table(probability_matrix(cgh), paste0(outdir_prefex,"_probability_matrix.txt"), sep='\t', quote=FALSE, col.names=NA)
write.table(copynumber(cgh), paste0(outdir_prefex,"_copynumber.txt"), sep='\t', quote=FALSE, col.names=NA)
write.table(calls(cgh), paste0(outdir_prefex,"_calls.txt"), sep='\t', quote=FALSE, col.names=NA)
write.table(segmented(cgh), paste0(outdir_prefex,"_segmented.txt"), sep='\t', quote=FALSE, col.names=NA)

sink()

