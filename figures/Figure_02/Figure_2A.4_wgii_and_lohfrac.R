################################################################################################################################################
##                                                                                                                      
##  SCRIPT TO GET LOG AND WGII
##                                                                                                                      
##  
##  Author: Hang Xu                                                                                                                
##           
##                                                                                                                      
################################################################################################################################################



library(dplyr)
get.loh.frac <- function(x){

	uni.chr <- unique(x$chromosome)
	a <- 0
	b <- 0

	for (chr in uni.chr){

		x.chr <- subset(x, chromosome==chr)
		x.loh.chr <- subset(x.chr, loh=="TRUE")

		b <- b + sum(as.numeric(x.chr[,"end"]) - as.numeric(x.chr[,"start"]))

		if (!is.null(dim(x.loh.chr))){
				a  <- a + sum(as.numeric(x.loh.chr[,"end"]) - as.numeric(x.loh.chr[,"start"]))
		}else{
			a <- a + sum(as.numeric(x.loh.chr[3]) - as.numeric(x.loh.chr[2]))
		}
	}
	return(round(a/b,4))
}

get.wgii <- function(cns){

	uni.chr <- unique(cns$chromosome)
	gii.weight <- 0

	for (chr in uni.chr){

		cns.chr <- subset(cns, chromosome==chr)
		cns.cn.chr <- subset(cns.chr, status!="NEU")

		b <- sum(as.numeric(cns.chr[,"end"]) - as.numeric(cns.chr[,"start"]))

		if (!is.null(dim(cns.cn.chr))){
			gii.chr <- sum(as.numeric(cns.cn.chr[,"end"]) - as.numeric(cns.cn.chr[,"start"]))/b
		}else{
			gii.chr <- sum(as.numeric(cns.cn.chr[3]) - as.numeric(cns.cn.chr[2]))/b
		}
		gii.weight <- gii.weight + gii.chr
	}
	return(round(gii.weight/length(uni.chr),4))
}

get.cnv <- function(cnv.fn, pp.fn){
	cnv <- read.delim(cnv.fn, stringsAsFactors=FALSE)
	pp <- read.delim(pp.fn, stringsAsFactors=FALSE)

	ploidy <- as.numeric(pp[1,"ploidy"])

	# remove sex chrom

	cnv <- subset(cnv, chromosome %in% paste0("chr",1:22))

	cnv$loh <- apply(cnv, 1, function(x) ifelse(round(as.numeric(x["minorAlleleCopyNumber"]))==0, "TRUE", "FALSE") )

	cnv$status <- NA

	for (j in 1:nrow(cnv)){

		if (as.numeric(cnv[j,"copyNumber"]) > 3*ploidy){
			cnv[j,"status"] <- "AMP"
		}else if (round(as.numeric(cnv[j,"minorAlleleCopyNumber"])) == 0 && round(as.numeric(cnv[j, "majorAlleleCopyNumber"])) == 0){  
			# this is equivalent to copyNumber < 0.5 in the data
			cnv[j,"status"] <- "HD"
		}else if (as.numeric(cnv[j,"copyNumber"]) - ploidy > 0.5){
			cnv[j,"status"] <- "GAIN"
		}else if (ploidy - as.numeric(cnv[j,"copyNumber"]) > 0.5){
			cnv[j,"status"] <- "LOSS"
		}else{
			cnv[j,"status"] <- "NEU"
		}
	}

	return(cnv)
}

frac.df <- c()

purple.dir <- "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_2/PURPLEV3"
setwd(file.path(purple.dir))
samples <- list.files(purple.dir)

`%notin%` <- Negate(`%in%`)

used <- list("wgii_and_LOHfrac.tsv")
for (sample.id in samples){
  sample_simple <- strsplit(sample.id, split='.purple.', fixed=TRUE)[[1]][1]
  print (sample_simple)
  if (sample_simple %notin% used){
  cnv.fn <- file.path(purple.dir, paste0(sample_simple, ".purple.cnv.somatic.tsv"))
  pp.fn <- file.path(purple.dir, paste0(sample_simple, ".purple.purity.tsv"))
  
  cnv <- get.cnv(cnv.fn=cnv.fn, pp.fn=pp.fn)
  loh.frac <- get.loh.frac(cnv)
  wgii <- get.wgii(cnv)
  frac.df <- rbind(frac.df, c(sample_simple, loh.frac, wgii))
  used <- c(used,sample_simple)
  }
}


colnames(frac.df) <- c("sample", "LOH.frac", "wGII")
write.table(frac.df, file="wgii_and_LOHfrac.tsv", sep="\t", quote=FALSE, row.names=FALSE)
