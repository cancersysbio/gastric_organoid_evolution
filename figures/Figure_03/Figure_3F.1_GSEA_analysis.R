#############################################################################################################################
##                                                                                                                      
##  GSEA ANALYSIS
##  FOR FIGURES 3F, 6G, ED7G, ED8F, ED9F
##                                                                                                                      
##  Date: 27 DECEMBER 2021                                                                                                                   
##  
##  Author: Kasper Karlsson
##
##                                                                                                                      
############################################################################################################################


### USE R version 3.6!

version

library(HTSanalyzeR2)
library(org.Hs.eg.db)
library(KEGGREST)
library(igraph)
library(GO.db)

### SET WORKING DIRECTORY

setwd("/Users/kasperkarlsson/_Stanford/Papers/Evolution_paper/Script_and_figures_final/Figure_scripts_commented/Figure_3/GSEA_and_Venn")



### FUNCTION TO GET GSEA ENRICHMENT FOR SEURAT OUTPUT FROM FIND MARKERS

geneSetAnalysis <- function(dfile) {
  dfile_subs <- dfile
  phenotype <- as.vector(dfile_subs$avg_log2FC) ### USE LOG FOLD CHANGE AS PHENOTYPE VECTOR
  names(phenotype) <- rownames(dfile_subs)

  ## specify the gene sets type you want to analyze
  MSig_H <- MSigDBGeneSets(species = "Hs", collection = "H", subcategory = NULL) # HALLMARK GENE SETS
  ListGSC <- list(MSig_H=MSig_H)
  
  ## iniate a *GSCA* object
  gsca <- GSCA(listOfGeneSetCollections=ListGSC, 
               geneList=phenotype)
  
  ## preprocess
  gsca1 <- preprocess(gsca, species="Hs", initialIDs="SYMBOL",
                      keepMultipleMappings=TRUE, duplicateRemoverMethod="max",
                      orderAbsValue=FALSE)
  
  ## analysis
  if (requireNamespace("doParallel", quietly=TRUE)) {
    doParallel::registerDoParallel(cores=4)
  }  
  
  gsca2 <- analyze(gsca1, 
                       para=list(pValueCutoff=0.05, pAdjustMethod="BH",
                                 nPermutations=10000, minGeneSetSize=1,
                                 exponent=1), 
                       doGSOA = FALSE)
  return(getResult(gsca2)$GSEA.results$MSig_H)
      message("Successfully ran GSEA.")
}

### READ IN FILE NAMES
realnames <- list.files(".", pattern="markers*", full.names=FALSE) # FILE NAMES WITHOUT PATH
filenames <- list.files(".", pattern="markers*", full.names=TRUE)  # FILE NAMES WITH PATH

ldf <- lapply(filenames,read.csv,sep="\t",header=TRUE,row.names=1) # READ IN ALL FILES
res <- lapply(ldf, try(geneSetAnalysis)) # GET GENE SET ENRICHMENT FOR EACH DEG FILE

realnames2 <- substr(realnames,1,nchar(realnames)-4) # GET SAMPLE NAME

### PRINT TO FILE
for (i in 1:length(res)){
  write.table(res[[i]],paste("GSEA/",realnames2[[i]],"_Hallmark_GSEA.txt",sep=""),sep="\t",quote=FALSE,row.names = TRUE,col.names=TRUE)
}
