################################################################################################################################################
# Population growth curves for non-ECB organoid cultures (all 3 donors), running from Early to Mid and Early to Late
#
# Estimate cumulative population size over time from cell counts and fit curve
# Examine growth derivatives for entire culture (non-ECB)
# 
#
# Eran Kotler
# Last updated 2021-05-25
# R version 4.0.3
#
################################################################################################################################################


library(tidyverse)
library(data.table)
require(gridExtra)
library(parallel)

source("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/sequtils/Supporting_programs/popGrowthHelper.R")



base_data_d <- "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_3/Figure_3B_Growth_trajectory_nonECB"
growth.data.d <- file.path(base_data_d, "Final_Input")
logbook.noECB.EarlyMid.f <- file.path(growth.data.d, 'Organoid_logbook_Mid_Trajectory.csv')
logbook.noECB.EarlyLate.f <- file.path(growth.data.d, 'Organoid_logbook_Late_Trajectory_Adjusted.csv') # From table 2, D1C1 and D2C1 thaw removed due to too few passages
results.d <- file.path(base_data_d, "Final_Output")



# Outputs:
deriv.by.samp.noECB.EarlyMid.f <- file.path(results.d, "deriv_by_samp_EarlyToMid.csv")
deriv.by.samp.noECB.EarlyLate.f <- file.path(results.d, "deriv_by_samp_EarlyToLate.csv")

EarlyMid.passage.csv.pfx <- file.path(results.d, "passage_num_conversion_EarlyMid")
EarlyLate.passage.csv.pfx <- file.path(results.d, "passage_num_conversion_EarlyLate")

EarlyMid.growth_curves_f <- file.path(results.d, "growth_curves_EarlyToMid.pdf")
EarlyLate.growth_curves_f <- file.path(results.d, "growth_curves_EarlyToLate.pdf")

# Change passage number to be consecutive (i.e. skip  missing passage information) and relative (start from 1)
adjust_passage <- function(df){
  passages <- df$Passage
  first_pass <- 1 
  num_pass <- length(passages)
  new_pass <- seq(as.numeric(first_pass), num_pass)
  return(as.numeric(new_pass))
}

organize_sample <- function(df_all, sample){
  # extract sample
  df <- df_all %>% filter(Sample==sample)
  df$Passage <- as.numeric(df$Passage)
  df <- df[order(df$Passage),] # Sort by passage
  df$Passage_corrected <- adjust_passage(df) # correct passage numbers (fix gaps)
  df$NumDays  <- df$Passage_corrected * 14 # Defined number of days by passage number (assuming 14 days between passages)
  return(df) 
}

popAtTime <- function(d, num.days, method='loess'){
  # calc growth curve, and fit approximation ('loess' / 'splines') to get approx. pop. size at time num.days
  d <- as.data.table(d)
  d <- d[!is.na(Total.viable.cells),]
  setorder(d, NumDays)
  d[, FCgrowth := Total.viable.cells / shift(Cells.seeded, 1L, type="lag")] # Growth (FC in size) from previous tp
  
  # Calc cumulative population size 
  d[1, "Cum.pop.size"] = d[1,"Total.viable.cells"]
  for (t in seq(2,nrow(d))){
    prev.pop <- d[t-1, "Cum.pop.size"]
    curr.FCgrowth <- d[t, "FCgrowth"]
    d[t, "Cum.pop.size"] <- curr.FCgrowth * prev.pop
  }
  d$Log.pop.size <- log10(d$Cum.pop.size+1) # Add 1 for numeric stability in log
  
  if (method=='loess'){ # fit LOESS and get value
    pop <- loess.approximation(d$NumDays, d$Log.pop.size, loess.span=2, x.to.pred = num.days,
                               show.plot=F, return_deriv=F)        
  }else {if (method=="splines"){ # fit splines and get value
    pop <- splines.approximation(d$NumDays, d$Log.pop.size, deg.fr=5, x.to.pred = num.days, 
                                 show.plot=F,return_deriv=F) 
  }}
}

plot_heatmap <- function(merged, title_str=""){
  # Plot heatmap (jsut for sanity check, not for publication)
  options(repr.plot.width = 8, repr.plot.height = 4)
  hm <- merged %>% 
    as.matrix %>% 
    #         sweep(., 2, .[1,]) %>% # Subtract first row to show change from T1
    as.data.frame() %>% 
    pivot_longer(-c(Passage), names_to = "samples", values_to = "Deriv") %>% # melt to long format
    ggplot(aes(x=Passage, y=samples , fill=Deriv)) + 
    geom_raster() +
    ggtitle(title_str)
  
  return(hm)
}

save_passage_numbers_csv <- function(samp_dfs, passage_csv_f_pfx){
  # Save csv file contating original passage number and corresponding newly-defined passage number for each sample
  samp_names <- names(samp_dfs)
  sapply(samp_names, function(samp_name){
    df <- samp_dfs[[samp_name]]
    df <- df[, c("Passage", "Passage_corrected")]
    colnames(df) <- c("Original passage number", "Adjusted passage number")
    out_csv_f <- paste0(passage_csv_f_pfx, "_", samp_name, ".csv")
    print (out_csv_f)
    write_csv(df, out_csv_f)           
    print(paste("Passage conversion table saved to", out_csv_f))
  })
}


make_sample_by_TP_derivs_table <- function(logbook_f, out_f, loess.span=2, plot_curves=FALSE, plots_pdf_f=NA, passage_csv_f_pfx=NA){
  # Main analysis function - load lab logbook counts data, plots growth curves and derivs for each sample
  # and saves a csv file with derivs for each sample across passages.
  # logbook_f: File containing count and seeding info for all samples in experiment
  # out_f: csv file to save table (for plotting in python/other)
  # loess.span: degree of smoothness of LOESS regression (larger -> more smooth)
  # plot_curves: If TRUE - plot growth curves and deirvatives for each sample
  # plots_pdf_f: path to pdf file onto which growth curves will be saved. If NA (and plot_curves==TRUE), curves are plotted to screen instead.
  # passage_csv_f_pfx: prefix (path basename) for csv files with conversion tables of Original passange number and newly created passage numbers (see note below).
  #                    A single file for each sample is saved (used if not NA). 
  #
  # Notes:
  #    - Time for growth curves & derivs is calculated by passage number X 14 (assuming two weeks between passages)
  #    - Passage numbers used are relative: renumbering passages in each sample from 1 -> end, in increments of 1, even where there 
  #      are missing info for a few passages in between. This slightly biases the rates of change of the derivatives at some points
  #      but was found more robust than imputation.
  
  ## Load & organize data 
  all <- read.csv(logbook_f, header=TRUE)
  all$Date2 <-  as.Date(all$Date,"%m/%d/%Y")
  startdate <- as.Date("2/14/17","%m/%d/%Y") 
  all$Sample <- paste0("D",all$Donor,"C",all$Culture)
  all <- all[all$Total.viable.cells!="Na",] # Remove rows with not viable cell count
  all <- all[all$Sample!="NA__",] # Remove empty rows 
  all$Total.viable.cells <- as.numeric(as.character(all$Total.viable.cells))
  all$Cells.seeded <- as.numeric(as.character(all$Cells.seeded))
  
  # Organize sample dataframe (passages, num of days, etc)
  samp_list <- unique(all$Sample) # Names of samples
  samp_dfs = lapply(samp_list, function(samp) {organize_sample(all, samp)})
  names(samp_dfs) <- samp_list
  print (samp_dfs)
  
  
  if (!is.na(passage_csv_f_pfx)){  # Save csv file with original and adjusted passage numbers
    save_passage_numbers_csv(samp_dfs, passage_csv_f_pfx)
  }
  
  # plot growth curves for each sample
  if (plot_curves==TRUE){
    options(warn=-1)
    if (!is.na(plots_pdf_f)){
      pdf(file = plots_pdf_f)
      print(paste("Growth curve plots saved to:", plots_pdf_f))
    }
    curves <- lapply(names(samp_dfs), function(samp) {growthCurve(samp_dfs[[samp]], title_pfx = samp, fit.type="loess", loess.span = loess.span)})    
    dev.off()
  }
  
  res <- mclapply(samp_dfs, function(d){ 
    d$Log.pop.size <- popAtTime(d, num.days=d["NumDays"], method='loess') # Get estimated log population size over time course
    samp_df <- tibble(Days=d$NumDays, Passage=d$Passage_corrected, LogNumCells=as.numeric(d$Log.pop.size))
    samp_df <- samp_df[is.finite(samp_df$LogNumCells),] 
    samp_df$Deriv <- sapply(samp_df$Days, function(t) {
      p <- loess.approximation(samp_df$Days, samp_df$LogNumCells, loess.span=loess.span, 
                               show.plot=F, x.to.pred=t, return_deriv=T) 
      p <- ifelse(!is.na(p), yes=p, 
                  no=loess.approximation(samp_df$Days, samp_df$LogNumCells, loess.span=loess.span, 
                                         show.plot=F, x.to.pred=t+1, return_deriv=T)) # if deriv is NA, get next day
    })
    
    return(samp_df[,c("Passage", "Deriv")]) 
  })
  
  # Merge all samples & organize
  merged <- Reduce(function(x, y) merge(x, y, by="Passage", all=TRUE), res)
  colnames(merged) <- c("Passage", samp_list)
  rownames(merged) <- merged$Passage
  
  # Save output table (csv) of derivative per sample:
  write_csv(merged, out_f)
  print(paste("Done, saved csv file to", out_f))
  
  return(merged)          
}

span <- 2 



early_mid <- make_sample_by_TP_derivs_table(logbook_f = logbook.noECB.EarlyMid.f, out_f = deriv.by.samp.noECB.EarlyMid.f, loess.span = span, plot_curves = TRUE,
                                             plots_pdf_f = EarlyMid.growth_curves_f,
                                             passage_csv_f_pfx = EarlyMid.passage.csv.pfx)

early_late <- make_sample_by_TP_derivs_table(logbook_f = logbook.noECB.EarlyLate.f, out_f = deriv.by.samp.noECB.EarlyLate.f, loess.span = span, plot_curves = TRUE,
                                             plots_pdf_f = EarlyLate.growth_curves_f,
                                             passage_csv_f_pfx = EarlyLate.passage.csv.pfx)

plot_heatmap(early_late, title_str = paste("Early-Late Derivs, LOESS span =",span))

