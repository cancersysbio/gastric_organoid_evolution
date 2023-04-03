################################################################################################################################################
##                                                                                                                      
##  GET NR OF KATAEGIS EVENTS USING MAFTOOLS
##                                                                                                                      
##  
##  Author: Kasper Karlsson                                                                                                                  
##           
##                                                                                                                      
################################################################################################################################################

library(maftools)

gorg_project_path <- "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/Revision/WGS_multiallel/maf_files_correct/"
setwd(file.path(gorg_project_path))

### Get nr kataegis events

D2C2_115d <- read.maf("D2C2_115d_Mutect2_multiallele_pass.vcf.maf")
D2C2_190d <- read.maf("D2C2_190d_Mutect2_multiallele_pass.vcf.maf")
D2C2_260d <- read.maf("D2C2_260d_Mutect2_multiallele_pass.vcf.maf")
D2C2_428d <- read.maf("D2C2_428d_Mutect2_multiallele_pass.vcf.maf")
D2C2_729d_LATE <- read.maf("D2C2_729d_LATE_Mutect2_multiallele_pass.vcf.maf")
D2C3_296d_MID <- read.maf("D2C3_296d_MID_Mutect2_multiallele_pass.vcf.maf")
D2C3_743d_LATE <- read.maf("D2C3_743d_LATE_Mutect2_multiallele_pass.vcf.maf")
D3C1_115d <- read.maf("D3C1_115d_Mutect2_multiallele_pass.vcf.maf")
D3C1_173d <- read.maf("D3C1_173d_Mutect2_multiallele_pass.vcf.maf")
D3C1_264d <- read.maf("D3C1_264d_Mutect2_multiallele_pass.vcf.maf")
D3C1_404d <- read.maf("D3C1_404d_Mutect2_multiallele_pass.vcf.maf")
D3C1_718d <- read.maf("D3C1_718d_LATE_Mutect2_multiallele_pass.vcf.maf")
D3C2_296d <- read.maf("D3C2_296d_MID_Mutect2_multiallele_pass.vcf.maf")
D3C2_750d <- read.maf("D3C2_750d_LATE_Mutect2_multiallele_pass.vcf.maf")
D3C3_296d <- read.maf("D3C3_296d_MID_Mutect2_multiallele_pass.vcf.maf")
D3C3_756d <- read.maf("D3C3_756d_LATE_Mutect2_multiallele_pass.vcf.maf")

rainfallPlot(maf=D2C2_115d, detectChangePoints = TRUE) # 0
rainfallPlot(maf=D2C2_190d, detectChangePoints = TRUE) # 1
rainfallPlot(maf=D2C2_260d, detectChangePoints = TRUE) # 1
rainfallPlot(maf=D2C2_428d, detectChangePoints = TRUE) # 1
rainfallPlot(maf=D2C2_729d_LATE, detectChangePoints = TRUE) # 1
rainfallPlot(maf=D2C3_296d_MID, detectChangePoints = TRUE) # 3
rainfallPlot(maf=D2C3_743d_LATE, detectChangePoints = TRUE) #3
rainfallPlot(maf=D3C1_115d, detectChangePoints = TRUE) # 1
rainfallPlot(maf=D3C1_173d, detectChangePoints = TRUE) # 0
rainfallPlot(maf=D3C1_264d, detectChangePoints = TRUE) # 1
rainfallPlot(maf=D3C1_404d, detectChangePoints = TRUE) # 2
rainfallPlot(maf=D3C1_718d, detectChangePoints = TRUE) # 2
rainfallPlot(maf=D3C2_296d, detectChangePoints = TRUE) # 2
rainfallPlot(maf=D3C2_750d, detectChangePoints = TRUE) # 1
rainfallPlot(maf=D3C3_296d, detectChangePoints = TRUE) # 1
rainfallPlot(maf=D3C3_756d, detectChangePoints = TRUE) # 1

