
# Environement
rm(list=ls())
graphics.off()
set.seed(1)


library(kableExtra)
library(ggplot2)
library(cowplot)
library(gGnome)
library(GenomicRanges)

jab.dir <- file.path("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_2/Jabba")
setwd(jab.dir)

sample = "D2C2_115d"
sample = "D2C2_190d"
sample = "D2C2_260d"
sample = "D2C2_428d"
sample = "D2C2_LATE"
sample = "D2C3_LATE"
sample = "D2C3_MID"
sample = "D3C1_115d"
sample = "D3C1_173d"
sample = "D3C1_264d"
sample = "D3C1_404d"
sample = "D3C1_LATE"
sample = "D3C2_LATE"
sample = "D3C2_MID"
sample = "D3C3_LATE"
sample = "D3C3_MID"
sample = "D3C1_LATE_server"
sample = "D3C1_LATE_bingxin"

gg <- readRDS(paste0("jabba.simple.gg_",sample,".rds"))
pdf(paste0("plot_3.6/",sample,"3:6.0e7-6.1e7.pdf"), width=20, height=10) # ALL D3 SAMPLES
plot(gg$gt,"3:6.0e7-6.1e7")
dev.off()

pdf(paste0("plot_3.8/",sample,"3:8.6e7-8.8e7.pdf"), width=20, height=10) # ALL D2 SAMPLES
plot(gg$gt,"3:8.6e7-8.8e7")
dev.off()
