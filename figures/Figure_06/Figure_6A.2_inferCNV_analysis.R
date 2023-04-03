################################################################################################################################################
##                                                                                                                      
##  Run infercnv                                                                             
##                                                                                                                      
##  Date: 11 March 2020                                                                                                                    
##  
##  Author: Moritz Przybilla                                                                                                                    
##                                                                                                                      
##                                                                                                                      
################################################################################################################################################

# clear workspace beforehand
rm(list = ls())

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("Matrix", "infercnv")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

# get arguments following the script
args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=1){
  message("\nError!\nusage: Rscript infercnv_analysis.R /path/to/patient_folder\n")
  quit()
}

# get working directory
w.dir <- args[1] # w.dir <- "/labs/ccurtis2/mjprzy/infercnv_gastric"
setwd(w.dir)

# get prepared input files 
cell.annotation.name <- paste0("cellAnnotations_", basename(w.dir), ".txt")
counts.sparse.matrix.name <- paste0("sc.10x.counts_",basename(w.dir),".RData")

# define output directory and create it
o.dir <- file.path(w.dir, "results")
dir.create(o.dir)

# read in prepared input matrix with positive and negative (reference) cells
load(counts.sparse.matrix.name)

# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix= mat.matrix.sparse,
                                    annotations_file= cell.annotation.name,
                                    delim="\t",
                                    gene_order_file= paste0(w.dir,"/gene_ordering_file.txt"),
                                    ref_group_names=c("normal"))

# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj, 
                             output_format="pdf",
                             num_threads = 16,
                             cutoff= 0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir= o.dir,  # dir is auto-created for storing outputs
                             window_length = 201, # 51 is way to small
                             cluster_by_groups=TRUE, # cluster
                             HMM=TRUE, # turn on to auto-run the HMM prediction of CNV levels
                             HMM_transition_prob=1e-6,
                             HMM_report_by = c("subcluster"), #,"cell"),
                             analysis_mode = c('subclusters'), #, 'cells'),
                             denoise=T,
                             sd_amplifier=1.35)  # sets midpoint for logistic
                      

# apply median filtering adding on to the one before
infercnv_obj_medianfiltered = infercnv::apply_median_filtering(infercnv_obj)

infercnv::plot_cnv(infercnv_obj_medianfiltered,
                   out_dir= o.dir,
                   output_format="pdf",
                   output_filename='infercnv.median_filtered',
                   x.range="auto",
                   x.center=1,
                   title = "infercnv",
                   color_safe_pal = FALSE)

