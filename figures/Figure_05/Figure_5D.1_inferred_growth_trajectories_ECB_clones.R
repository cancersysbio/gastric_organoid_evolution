################################################################################################################################################
# RG-specific "absolute fitness" estimation by growth curve analysis on a specific culture.
#
# Main steps:
# 1) loads and organizes lab logbook info (cell counts, splits, etc.) and RG barcode DNA-seq results (relative fractions fo each RG)
# 2) Creates a bulk growth curve for the entire culture 
# 3) Calculates RG-specific growth curves using their relative fraction in the popoulation at each time point
# 4) Fits a line through absolute abundances of the RGs and estimates the derivative for each RG at each time point (using LOESS or splines)
#
# Outputs are:
# - Sample (bulk) growth curve) - only displayed)
# - PDF files with RG-specific growth curves for each Rep (top n RGs are plotted)
# - CSV file with absolute population sizes
# - A CSV file for each rep with estimated absolute fitness at each time point (derivatives)
# - A merged CSV file with fitnesses of all 3 replicates together
#
# Run entire workflow using the function doGrowthCurveAnalysis()
#
# Eran Kotler
# Last updated 2021-05-25
# R version 4.0.3
#
################################################################################################################################################



load_packages <- function(){
	pkgs_needed = c("tidyverse","dplyr","ggplot2","data.table","parallel")

	if (!requireNamespace("BiocManager", quietly = TRUE))
	    install.packages("BiocManager")
	BiocManager::install(setdiff(pkgs_needed, installed.packages()))

	suppressPackageStartupMessages({
	    sapply(pkgs_needed, require, character.only = TRUE)
	    })

	source("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/sequtils/Supporting_programs/popGrowthHelper.R") # Helper functions for growth curves
}


set_globals <- function(samp_time, clone){
	# Set global variables for analysis

	assign("SAMP_TIME", samp_time, envir = .GlobalEnv)
	assign("CLONE", clone, envir = .GlobalEnv)

	if (clone=="D2C1"){
		assign("PATIENT", "D2", envir = .GlobalEnv)
	}
	if (clone=="D2C2"){
		assign("PATIENT", "D2", envir = .GlobalEnv)
	}
	if (clone=="D2C3"){
		assign("PATIENT", "D2", envir = .GlobalEnv) 
	}

	if (samp_time=="Early" & clone=="D3C2"){
		assign("PATIENT", "D3", envir = .GlobalEnv) 
		assign("CLONE", "D3C2", envir = .GlobalEnv)

	}

	print("Global vars set")
}

get_paths <- function(samp_time, clone){
	# Paths to input and output files
	base_data_d <- "/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_5/Growth_trajectory_ECB/"
	general_paths <- list( # Common to all samples
		base_data_d=base_data_d,
		growth_data_d=file.path(base_data_d, "Input/Lab_logbooks"),
		barcode_data_d=file.path(base_data_d, "Input/norm_RGC_matrices"),
		results_d="Output"
		)

	# Specific paths (for individual experiments/clones etc.
	if (samp_time=="Early"){
		logbook_ECB_f <- file.path(general_paths$growth_data_d, "Organoid logbook ECB Simple.csv")

		if (clone=="D2C1"){
			RGC_mat_f <- file.path(general_paths$barcode_data_d, "D2C1_RGC_matrix_norm.txt") }
		else{if (clone=="D2C2"){
			RGC_mat_f <- file.path(general_paths$barcode_data_d, "D2C2_RGC_matrix_norm.txt") } 
		else{if (clone=="D2C3"){
			RGC_mat_f <- file.path(general_paths$barcode_data_d, "D2C3_RGC_matrix_norm.txt") } 
		else{if (clone=="D3C2"){
			RGC_mat_f <- file.path(general_paths$barcode_data_d, "D3C2_RGC_matrix_norm.txt") } 
			}}}
		}

	specific_paths <- list(
		logbook_ECB_f=logbook_ECB_f, 
		RGC_mat_f=RGC_mat_f
	)

	# Outputs:
	RG_abs_pop_sizes_s <- sprintf("%s_%s_RG_abs_pop_sizes.csv", clone, samp_time) # For storing absolute population sizes
	RG_abs_fitness_base_s <- sprintf("%s_%s_RG_abs_fitness", clone, samp_time) # Basename (path prefix) for storing estimated absolute fitness

	outputs <- list(
		RG_abs_fitness_base=file.path(general_paths$results_d, RG_abs_fitness_base_s),
		RG_abs_pop_sizes_f=file.path(general_paths$results_d, RG_abs_pop_sizes_s)
		)

	paths <- c(general_paths, specific_paths, outputs)
	# print(paths)
	return(paths)
}

get_logbook_data <- function(paths, samp_time){
	# Load lab notebook information (cell counts, seedings, dilutions, etc.)
	
	if (samp_time=="Early"){skip.rows <- 10}else{if 
		(samp_time=="Late"){skip.rows <- 0}	}

	lb <- read.csv(paths$logbook_ECB_f, header=TRUE, skip = skip.rows) # skipping first rows (describe # sorted cells)

#	lb <- lb  %>% filter(!is.na(Cells.seeded))
	lb <- subset(lb,Cells.seeded != "Na")
	lb <- subset(lb,Total.viable.cells != "Na")
	lb$Total.viable.cells <- as.numeric(gsub(",","",lb$Total.viable.cells))
	lb$Cells.per.well <- as.numeric(gsub(",","",lb$Cells.per.well))
	lb$Cells.seeded <- as.numeric(gsub(",","",lb$Cells.seeded))
	startdate <- as.Date("2/14/17","%m/%d/%Y") 
	lb$Date2 <- as.Date(lb$Date,"%m/%d/%Y")
	lb$NumDays <- as.numeric(difftime(lb$Date2, startdate ,units="days"))
	lb$Sample <- as.factor(paste(lb$DC, lb$Replicate, sep="_"))
	lb$Dilution.factor <- lb$Total.viable.cells / lb$Cells.seeded # dilution factor in current seeding
	return(lb)
}

# get estimated total pop size at each TP
# Add log10(total-pop-size) for each clone at each time point:
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
    d$Log.pop.size <- log10(d$Cum.pop.size+1) # +1 for numeric stability in log
    if (method=='loess'){ # fit LOESS and get value
        pop <- loess.approximation(d$NumDays, d$Log.pop.size, loess.span=2, x.to.pred = num.days,
                                   show.plot=F, return_deriv=F)        
    }else {if (method=="splines"){ # fit splines and get value
        pop <- splines.approximation(d$NumDays, d$Log.pop.size, deg.fr=5, x.to.pred = num.days, 
                                     show.plot=F,return_deriv=F) 
    }}
}

get.pop <- function(data, samp, method='loess'){
    d <- data[data$Sample==samp["Sample"],]
    t <- samp["NumDays"]
    popAtTime(d, num.days=as.numeric(t), method=method)
}

load_rg_info <- function(RGC_f){
	# Get barcode-seq data from RGC file
	rg <- fread(RGC_f)
	colnames(rg)[1] <- "RG"
	rownames(rg) <- rg$RG
	return(rg)
}

make_samp_metadata <- function(rg){
	# create samp info (metadata) for RG column names (for simple merging later)
	samp.met <- data.frame(row.names = colnames(rg)[-1], rg.col.name = colnames(rg)[-1])

	# pat.str.pos <- ifelse(samp_time=="Early", 1, 2) # position of Patinet in string
	# samp.met$patient <- as.factor(sapply(as.character(samp.met$rg.col.name), function(x) strsplit(x, "_")[[1]][pat.str.pos]))
	# samp.met$clone <- as.factor(sapply(as.character(samp.met$rg.col.name), function(x) substr(x, start=regexpr("_C", x)+1, stop = regexpr("_C", x)+2)))
	samp.met$clone <- CLONE # Get clone from defined params rather than from rg column name
	samp.met$rep <- as.factor(sapply(as.character(samp.met$rg.col.name), function(x) substr(x, start=regexpr("_R", x)+1, 
	                                                                          stop = regexpr("R", x)+1)))
	samp.met$tp <- as.factor(sapply(as.character(samp.met$rg.col.name), function(x) tail(strsplit(x, "_")[[1]],1)))
	return(samp.met)
}

organize_lb <- function(lb, analyzed_clone){
	# Organize logbook data, add log-pop size estimation and organize columns for merging with RG data
  
  #print ("lb2",lb)
  no_parent <- filter(lb, !Replicate %in% c("PARENT1", "PARENT 2", "Parent","PARENT","PARENT 1")) # remove parent timepoint (common to all reps)
#  print ("no parent",no_parent)
	no_parent$Log.pop.size <- apply(no_parent, 1, function(x) get.pop(data=no_parent, samp=x, method='loess'))
#	print ("no parent2",no_parent)
	no_parent$samp.name2 <- apply(no_parent, 1, function (x){
	    paste(x["DC"], x["Replicate"], x["Time.point"], sep="_")})
  
	#print (no_parent)
	
	return (no_parent)
}

plot_culture_growth_curves <- function(lb, out_d){
	sapply(c("R1", "R2", "R3"), function(r) {
		d <-filter(lb, Replicate==r)
		fname <- sprintf("%s_%s_%s_culture_growth.pdf", CLONE, SAMP_TIME, r)
		out_f <- file.path(out_d, fname)
		title_str <- sprintf("%s %s %s growth", CLONE, SAMP_TIME, r)
		pdf(out_f, width=8, height=6, pointsize=0.1, )
		pl <- growthCurve(d, title_pfx=title_str, fit.type="loess", loess.span = 2) 
		dev.off()
		})
	}

merge_lb_and_meta <- function(lb, samp.met, clone){
	# Create a merged dataframe with lab logbook info and sample metadata (col names based) from RG info
	# lb should be a lab logbook (cell counts etc) df with no parent time point
	
  print ("lb")
  print (lb)
  print ("samp.met")
  print (samp.met)
  
  samp.met$samp.name2 <- apply(samp.met, 1, function (x) {
		paste(clone, x["rep"], x["tp"], sep="_") 
		})
  
  
  print ("samp.met2")
  print (samp.met)
  
	df <- merge(samp.met, lb, by="samp.name2", suffixes = c("","_"))
	cols.to.keep <- c('Sample', 'samp.name2', 'rg.col.name','clone','Organoid','Replicate','Viability',
	        'Total.viable.cells','Cells.per.well','Cells.seeded','Date','Time.point','Date2','NumDays','Dilution.factor',
	        'Log.pop.size') # patient
	df <- dplyr::select(df, all_of(cols.to.keep))

	return(df)
}

abs_fitness <- function(rg, df){
	# Create a dataframe with estimated "absolute fitness" for each RG at each time point by:
	# calculating: 
	# 		log(RG pop size at time TP) = log10(total pop at TP) + log10(relative fract of RG at TP)

	#Convert rg table (relative frequencies) to log10(rg) to comply with log10(pop size)
	rg.log10 <- cbind(RG=rg$RG, log10(rg[,!"RG"]))

	# For each RG calc: log(RG pop size at time TP) = log10(total pop at TP) + log10(relative fract of RG at TP)
	log.tot.pops <- df$Log.pop.size
	names(log.tot.pops) <- as.character(df$rg.col.name)
	log.rel.fracts <- dplyr::select(rg.log10, all_of(names(log.tot.pops)))
	abs.rg <- cbind(RG=rg.log10$RG, log.rel.fracts + log.tot.pops[col(log.rel.fracts)])
	return(abs.rg)
}

RGGrowthCurve <- function(abs.rg, RGnum, samp.names, df, title.pfx="", loess.span=2, deg.fr=3, x.to.pred=NA, method='loess'){
    # Create growth curve and derivs of line fit using LOESS (loess.span is importnat) or splines (deg.fr imporant).
    # if <x.to.pred> is specified - prints the value of the derivative at this x position

    d.rg <- abs.rg[RG==RGnum, ..samp.names] # extract data for single RG
    cnames <- names(d.rg) # get sample names
    days <- sapply(cnames, function (x) df[df$rg.col.name==x,"NumDays"]) # convert names to time
    a <- tibble(Days=days, LogNumCells=as.numeric(d.rg))
    a <- a[is.finite(a$LogNumCells),] # remove data points with no support (0 counts in seq data)
    if (method=='loess'){
        if (!is.na(x.to.pred)){
            s <- sprintf("Slope at %i: %f", x.to.pred,
                        loess.approximation(a$Days, a$LogNumCells, loess.span=loess.span, 
                                            show.plot=F, x.to.pred=x.to.pred, return_deriv=T))
           print(s)
        }

        p <- loess.approximation(a$Days, a$LogNumCells, loess.span=loess.span, show.plot=T, title_pfx=title.pfx)
    }else{if (method=='splines'){
        p <- splines.approximation(a$Days, a$LogNumCells, deg.fr=deg.fr, show.plot=T, title_pfx=title.pfx)
    }}
}

plot_top_RG_growthCurves_to_file <- function(abs.rg, samp.names, df, n_to_plot, out_f){
	# Save RG-specific growth-curve plots to pdf (<n_to_plot> top RGs are plotted to <out_f>)
	pdf(out_f, width=8, height=6, pointsize=0.1, )
	sapply(head(abs.rg$RG, n_to_plot), function(x) RGGrowthCurve(abs.rg, RGnum=x, samp.names=samp.names, df=df,
	                                                       title.pfx=x, method='loess',loess.span=2))
	dev.off()	
	print(sprintf("Growth curves saved to %s", out_f))
}

getSlopes <- function(abs.rg, RGnum, samp.names, df, method='loess', loess.span=2){
	# As a measure of absolute fitness - take derivative of log(population size) for an RG at a given time point
    d.rg <- abs.rg[RG==RGnum, ..samp.names] # extract data for single RG
    cnames <- names(d.rg) # get sample names
    days <- sapply(cnames, function (x) df[df$rg.col.name==x,"NumDays"]) # convert names to time
    a <- tibble(Days=days, LogNumCells=as.numeric(d.rg))
    a <- a[is.finite(a$LogNumCells),] # remove data points with no support (0 counts in seq data)

    if (dim(a)[1]<4){ # require min of 4 data points for interpolation
        return(rep(NA, times = length(days))) 
        } else {                    
            if (method=='loess'){
                slopes <- sapply(days, function(t) {
                    p <- loess.approximation(a$Days, a$LogNumCells, loess.span=loess.span, 
                                       show.plot=F, x.to.pred=t, return_deriv=T) 
                    p <- ifelse(!is.na(p), yes=p, 
                                no=loess.approximation(a$Days, a$LogNumCells, loess.span=loess.span, 
                                       show.plot=F, x.to.pred=t+1, return_deriv=T)) # if deriv is NA, get next day
                    })
                
                nms <- names(slopes)
                slopes <- as.numeric(slopes)
                names(slopes)=nms
                
                return(slopes)
                } 

            else{ if (method=='splines'){
                slopes <- sapply(days, function(t) {
                   splines.approximation(a$Days, a$LogNumCells, x.to.pred=t, deg.fr=3, show.plot=F, return_deriv=T) 
                })
            return(slopes)
    }}}
}

getRGslopesForSamps <- function(abs.rg, samples, times.df, loess.span=2){
	# Wrapper for getting slopes (derivs) for multiple samples (all RGs in each)
    non.exist.cols <- setdiff(samples, colnames(abs.rg))
    samples <- samples [! samples %in% non.exist.cols] # filter out irrlelvant samples (not in count data)
    rg.slopes <- sapply(abs.rg$RG, function(x) getSlopes(abs.rg, RGnum=x, samp.names=samples, df=times.df, 
                                                         method="loess", loess.span=loess.span))
    rg.slopes %>% t() %>% as.data.frame # organize
}

plot_top_RG_growthCurves_all_reps <- function(paths, clone, cl_reps, samp_time, dt, df, abs.rg, n_to_plot=10){
	# wrapper to create pdfs with RG-specific growth curves and derivs - a separate pdf for each replicate within the clone
	lapply(cl_reps, function(r){
		RG_rep_growth_f=file.path(paths$results_d, sprintf("RG_Cum_growth_loess_%s%s_%s.pdf", clone, r, samp_time))
		
		samp_names <- dt[rep==r & clone==clone, rg.col.name]
		samp_names <- intersect(samp_names, colnames(abs.rg)) # remove time points that exist in RG info but not in logbook (e.g. C1R2-T8)

		plot_top_RG_growthCurves_to_file(abs.rg=abs.rg, samp.names=samp_names, df=df, n_to_plot=10, out_f=RG_rep_growth_f)
		})
}

                        

# **** MAIN FUNCTION FOR GROWTH CURVE ANALYSIS ON A GIVEN CLONE *****
doGrowthCurveAnalysis <- function(samp_time, clone, cl_reps = c("R1","R2","R3"), loess_span=2){
	# <samp_time>: Experiment time from {"Early", "Late"}
	# <clone>: Clone to analyze (for Late can be one of {"D2C1", "D2C2" "D2C3"})
	# <loess_span>: span to use for LOESS line fitting/smoothing (default 2)
	# <cl_reps>: Names of replicates available within the specified clone (usually c("R1","R2","R3")) 
	# note: Current growth curves do not use "Parent" time points as these are not deignated to any replicate.
	
	set_globals(samp_time=samp_time, clone=clone) # Set global vars for given sample (PATIENT, CLONE, SAMP_TIME)
	print (sprintf("Starting analysis for %s %s", samp_time, CLONE))

	print (samp_time)
	paths <- get_paths(samp_time=samp_time, clone=CLONE) # Get file paths for inputs and outputs

	lb_data <- get_logbook_data(paths, samp_time) # Get lab logbook data (cell counts etc)
	print ("lb_data")
	print (lb_data)
	no_parent <- organize_lb(lb_data, clone)
	print ("no_parent")
	print (no_parent)

	plot_culture_growth_curves(no_parent, paths$results_d) # Plot growth dynamics and derivs for all reps

	rg <- load_rg_info(paths$RGC_mat_f) # Get barcode-seq data    
	samp.met <- make_samp_metadata(rg)
	head(no_parent)
	head(samp.met)
	df <- merge_lb_and_meta(no_parent, samp.met, CLONE) # Organized cell counts info
	abs.rg <- abs_fitness(rg, df) # Absolute rg abundance over time
	dt <- data.table(samp.met) # sample info data.table

	plot_top_RG_growthCurves_all_reps(paths, CLONE, cl_reps, samp_time, dt, df, abs.rg, n_to_plot=10)

	# Save absolute population sizes
	write_csv(abs.rg, paths$RG_abs_pop_sizes_f) 
	#print(sprintf("Absolute population sizes saved to %s", paths$RG_abs_pop_sizes_f))

	# Get RG slopes for all Clone1 reps in this dataset and save csv files
	cl_samp_ls <- lapply(cl_reps, function (rep_num) dt[rep==rep_num & clone==CLONE, rg.col.name]) # sep clone samples by rep

	res <- mclapply(cl_samp_ls, function(x) {
	    getRGslopesForSamps(abs.rg, samples=x, times.df=df, loess.span=loess_span)
	    })
	names(res) <- cl_reps

	lapply(cl_reps, function(r){
	    f <- paste0(paths$RG_abs_fitness_base, "_", r ,".csv")
	    write.csv(res[r], f, row.names=TRUE)
	    print(paste(r, "Saved to", f))
	})

	# Merge R1-R3 to a single dataframe and save
	unified = cbind(res[[1]], res[[2]], res[[3]])

	f <- paste0(paths$RG_abs_fitness_base, "_allReps.csv")
	write.csv(unified, f, row.names=TRUE)
	print(paste("Unified table (R1-R3) saved to", f))

	print ("DONE")
}

### ============================================================================================================================================

                         
                         
#################
# Running script:
#################
load_packages()

setwd(file.path("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_5/Growth_trajectory_ECB/"))
library("tidyverse")
library("dplyr")
library("ggplot2")
library("data.table")
library("parallel")
source("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/sequtils/Supporting_programs/popGrowthHelper.R")


doGrowthCurveAnalysis(samp_time="Early",clone="D2C1")
doGrowthCurveAnalysis(samp_time="Early",clone="D2C2")
doGrowthCurveAnalysis(samp_time="Early",clone="D2C3")
doGrowthCurveAnalysis(samp_time="Early",clone="D3C2")
