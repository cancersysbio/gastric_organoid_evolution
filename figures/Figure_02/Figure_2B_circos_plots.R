################################################################################################################################################
##                                                                                                                      
##  SCRIPT TO CREATE CIRCOS PLOTS
##                                                                                                                      
##  
##  Author: Hang Xu                                                                                                                
##           
##                                                                                                                      
################################################################################################################################################


library(circlize)
library(GenomicRanges)
library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(vcfR)


path = "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_2/"
setwd(file.path(path))

sample = "D3C1_115d"
snv.fn = paste0(path,"WGS_MUTECT2_multiallele_vaf/",sample,"_Mutect2_multiallele_pass.vcf")
cns.fn = paste0(path,"PURPLEV3/",sample,".cns")
sv.fn = paste0(path,"PURPLEV3/",sample,".purple.sv.vcf")
pdf.fn = paste0(path,"Circos/plots/",sample,"_circos.pdf")

# function to get SV type
get.svtype <- function(chr1, strand1, chr2, strand2){

    if (chr1==chr2){
        if (strand1=="+" && strand2=="-"){
            svtype <- "DEL"
        }else if (strand1=="-" && strand2=="+"){
            svtype <- "DUP"
        }else if (strand1=="+" && strand2=="+"){
            svtype <- "h2hINV"
        }else if (strand1=="-" && strand2=="-"){
            svtype <- "t2tINV"
        }else{
            stop(paste("Unknown svtype",strand1,strand2))
        }
    }else{
        svtype <- "TRA"
    }
    return(svtype)
}

# read in data

sv <- readVcf(sv.fn, "hg38")
sv.gr <- breakpointRanges(sv)
sv.bedpe <- breakpointgr2bedpe(sv.gr)

mut <- read.vcfR(snv.fn)
mut.tidy <- vcfR2tidy(mut, single_frame = TRUE)
mut.tidy.df <- subset(as.data.frame(mut.tidy$dat), Indiv!="" & !grepl("WT",Indiv) & !is.na(gt_AF))
mut.tidy.df$gt_AF <- as.numeric(mut.tidy.df$gt_AF)

cns <- read.delim(cns.fn, stringsAsFactors=F)
cns$width <- cns$end-cns$start
cns <- subset(cns, log2<=1.5 & log2>=(-1.5) & width>100000)

# get SV type
sv.bedpe$svtype <- apply(sv.bedpe, 1, function(x) get.svtype(x["chrom1"], x["strand1"], x["chrom2"], x["strand2"]))
sv.bedpe$svtype <- gsub("[h2h|t2t]", "", sv.bedpe$svtype)
sv.bedpe <- subset(sv.bedpe, chrom1 %in% paste0("chr",1:22) & chrom2 %in% paste0("chr",1:22))

svs <- as.data.frame(table(sv.bedpe$svtype))

write.csv(svs,paste0("Figure_2A/Structural_variations_",sample,".csv"),quote=FALSE) ### USED AS INPUT FOR STRUCTURAL VARIANTS IN FIGURE 2A

# format SVs into pairs
sv.pair.1 <- sv.bedpe[,c("chrom1","start1","end1")]
sv.pair.2 <- sv.bedpe[,c("chrom2","start2","end2")]

# SV colors
sv.palette <- c("#3FA7D6","#F97F6C","#FAC05E","#F17FAB")
names(sv.palette) <- c("DEL","DUP","INV","TRA")
sv.col <- sv.palette[sv.bedpe$svtype]

# plot circos
pdf(pdf.fn, width=6, height=6)

circos.initializeWithIdeogram(species = "hg38", plotType=c("labels","ideogram"), labels.cex=1.2, chromosome.index = paste0("chr",1:22))

circos.genomicTrack(mut.tidy.df[,c("CHROM","POS","POS","gt_AF")],
		numeric.column = 4,
		panel.fun = function(region, value, ...) {
    	circos.genomicPoints(region, value, pch=20, cex=0.2, col="#BB8538")
	},
	bg.border=NA
)
circos.genomicTrack(cns[,c("chromosome","start","end","log2")],
		numeric.column = 4,
		panel.fun = function(region, value, ...) {
    	circos.genomicLines(region, value, pch=20, lwd=2, col="#59CD90", type = "segment")
	},
	bg.border=NA
)
circos.genomicLink(sv.pair.1, sv.pair.2, col=sv.col, lwd=2)

dev.off()


# END #

