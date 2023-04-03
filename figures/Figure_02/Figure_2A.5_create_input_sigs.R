################################################################################################################################################
##                                                                                                                      
##  RUN HDP PROCESS FOR THE CLUSTERS DETERMINED BY NDP
##                                                                                                                      
##  Date: 19 AUGUST 2021                                                                                                                    
##  
##  Author: Moritz Przybilla                                                                                                                    
##           
##                                                                                                                      
################################################################################################################################################

# clear workspace beforehand
rm(list = ls())
options(stringsAsFactors = F)

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("hdp")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

#####################################################################################
# READ IN DATA
#####################################################################################

# read in the input
output.dir=commandArgs(T)[1]
n=as.numeric(commandArgs(T)[2])
mutations=read.table(paste0(output.dir, "/trinuc_mut_mat_full.txt"))
key_table=read.table(paste0(output.dir, "/key_table_full.txt"), header = T)

# remove everything below 100 mutations and count the total occurences per patient
# sample_remove=rownames(mutations)[rowSums(mutations)<100]
# mutations=mutations[!rownames(mutations)%in%sample_remove,]
# key_table=key_table[!key_table$Sample%in%sample_remove,]
freq=table(key_table$Patient)

#####################################################################################
# RUN THE HDP PROCESS
#####################################################################################

ppindex <- c(0, rep(1, nrow(mutations)))
cpindex <- c(1, rep(2, nrow(mutations)))

hdp_mut <- hdp_init(ppindex = ppindex, # index of parental node
                    cpindex = cpindex, # index of the CP to use
                    hh = rep(1, 96), # prior is uniform over 96 categories
                    alphaa = rep(1, length(unique(cpindex))), # shape hyperparameters for 2 CPs
                    alphab = rep(1, length(unique(cpindex))))  # rate hyperparameters for 2 CPs

hdp_mut <- hdp_setdata(hdp_mut, 
                       dpindex = 2:numdp(hdp_mut), # index of nodes to add data to
                       mutations) # input data (mutation counts, sample rows match up with specified dpindex)

hdp_activated <- dp_activate(hdp_mut, 1:numdp(hdp_mut), initcc=10, seed=n*200)

# Run posterior sampling chains
chain=hdp_posterior(hdp_activated,
                    burnin=20000,
                    n=100,
                    seed=n*1000,
                    space=200,
                    cpiter=3)

# save the chain results here
saveRDS(chain, paste0(output.dir, "/hdp_chain_",n,".Rdata"))
