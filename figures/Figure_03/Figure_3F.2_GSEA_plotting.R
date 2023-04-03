#############################################################################################################################
##                                                                                                                      
##  plot GSEA pathways for organoids evolution paper
##                                                                                                                      
##  
##  Author: Katie Houlahan
##
##                                                                                                                      
############################################################################################################################


### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(gridExtra)
library(tidyr)

setwd('/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_3/')

# set date
date <- Sys.Date()

### FIGURE 3J #####################################################################################
# read in gsea results 
files <- list.files(path = 'Figure_3F_GSEA', pattern = 'markers', full.name = TRUE)

gsea <- sapply(files, function(x) {
	tmp <- read.delim(x, as.is = TRUE)
	tmp[order(rownames(tmp)),]
	}, 
	simplify = FALSE)
plot_data_es <- do.call(cbind, lapply(gsea, '[[', 1))
plot_data_fdr <- do.call(cbind, lapply(gsea, '[[', 3))

colnames(plot_data_es) <- sapply(colnames(plot_data_es), function(x) paste(strsplit(x, '_')[[1]][3:4], collapse = '_'))
colnames(plot_data_fdr) <- sapply(colnames(plot_data_fdr), function(x) paste(strsplit(x, '_')[[1]][3:4], collapse = '_'))
rownames(plot_data_es) <- rownames(gsea[[1]])
rownames(plot_data_fdr) <- rownames(gsea[[1]])

# only keep a subset of pathways
pathways <- c('HALLMARK_TNFA_SIGNALING_VIA_NFKB','HALLMARK_HYPOXIA','HALLMARK_APOPTOSIS','HALLMARK_CHOLESTEROL_HOMEOSTASIS',
	'HALLMARK_MTORC1_SIGNALING','HALLMARK_P53_PATHWAY','HALLMARK_OXIDATIVE_PHOSPHORYLATION','HALLMARK_MYC_TARGETS_V1',
	'HALLMARK_E2F_TARGETS','HALLMARK_DNA_REPAIR','HALLMARK_G2M_CHECKPOINT')
plot_data_es <- plot_data_es[pathways,]
plot_data_fdr <- plot_data_fdr[pathways,]

samples <- c('D1C1_MID','D1C1_LATE','D1C2_MID','D1C2_LATE','D1C3_MID','D1C3_LATE',
	'D2C2_MID','D2C2_LATE','D2C3_MID','D2C3_LATE','D3C2_MID','D3C2_LATE')
plot_data_es <- plot_data_es[,samples]
plot_data_fdr <- plot_data_fdr[,samples]

# set yaxis lab 
yaxis.lab <- rownames(plot_data_es)
yaxis.lab <- gsub("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "TNFA", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_HYPOXIA", "HYPOXIA", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_APOPTOSIS", "APOPT.", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_CHOLESTEROL_HOMEOSTASIS", "CHOL HOM", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_MTORC1_SIGNALING", "MTORC1", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_P53_PATHWAY", "P53", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "OX PHOS", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_MYC_TARGETS_V1", "MYC V1", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_E2F_TARGETS", "E2F", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_DNA_REPAIR", "DNA REP", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_G2M_CHECKPOINT", "G2M", yaxis.lab)

# set spot size and colour functions
spot.size.function <- function(x) {abs(x)*3}
spot.colour.function <- function(x) {
    colours <- rep('white', length(x));
    colours[x > 0] <- 'firebrick3';
    colours[x < 0] <- 'dodgerblue3';
    colours[x == 0] <- 'transparent';
    return(colours);
    }

# set covariate to show clusters and cohort
time_cov <- rep(c("#18188B", "#6A2A58"), 6)
sample_cov <- rep(c("#BBDD34", "#00B050","#1E4620","#843C0C","#312113","#7C7C7C"), each = 2)

covariate <- list(
	rect = list(
		col = 'black',
		fill = time_cov,
		lwd = 1.5
		),
	rect = list(
		col = 'black',
		fill = sample_cov,
		lwd = 1.5
		)
	)

cov_grob <- covariates.grob(
	covariates = covariate,
	ord = length(time_cov):1,
	side = 'right'
	)

cov_legend <- list(
	legend = list(
		colours = c("#18188B", "#6A2A58"),
		labels = c('Mid','Late'),
		title = 'Timepoint'
		),
	legend = list(
		colours = c("#BBDD34", "#00B050","#1E4620","#843C0C","#312113","#7C7C7C"),
		labels = unique(gsub('_LATE|_MID', '', colnames(plot_data_es))),
		title = 'Culture'
		)
	)

legend_grob <- legend.grob(
	legends = cov_legend,
	label.cex = 1.8,
	title.cex = 1.8
	)

key_sizes <- seq(1, -1, -0.5);
key_cov <- draw.key(
	list(
			space = 'right',
			points = list(
				cex = spot.size.function(key_sizes),
				col = 'black',
				fill = spot.colour.function(key_sizes),
				pch = 21
				),
			text = list(
				lab = c('1.0','0.5','0.0','-0.5','-1.0'),
				cex = 1.5,
				adj = 1.5,
				fontface = 'bold'
				),
			title = expression(bold("Observed\nScore")),
			cex.title = 1.8,
			background = 'white',
			padding = 8
			)
	)


# create dotmap
create.dotmap(
		x = t(plot_data_es),
		bg.data = t(-log10(plot_data_fdr)),
		filename = 'Figure_3F_GSEA/Figure3F.pdf',
		xaxis.cex = 1.2,
		xaxis.rot = 45,
		yaxis.tck = 0,
		xaxis.tck = 0,
		spot.colour.function = spot.colour.function,
		spot.size.function = spot.size.function,
		xaxis.lab = paste('   ', yaxis.lab),
		yaxis.lab = rep('', ncol(plot_data_es)),
		colour.scheme = c('white','black'),
		na.spot.size = 0,
		bg.alpha = 1,
		colourkey = TRUE,
		colourkey.labels = c(
			expression(1),
			expression(10^-2),
			expression(10^-4),
			expression(10^-6)
			),
		at = c(0, seq(1.3, 6.1, 0.2)),
		colourkey.labels.at = seq(0, 6, 2),
		colourkey.cex = 1.5,
		axis.top = 1.5,
		key = NULL,
		legend = list(
			left = list(
				fun = cov_grob
				),
			inside = list(
				fun = legend_grob,
				x = -0.39,
				y = 1
				),
			inside = list(
				fun = key_cov,
				x = 1.05,
				y = 0.95
				),
			inside = list(
				fun = draw.key(list(text = list(lab = ''), title = expression(bold('FDR')), cex.title = 1.8)),
				x = -0.1,
				y = -0.29
				)
			),
		width = 11,
		height = 5,
		bottom.padding = 3,
		right.padding = 20,
		left.padding = 25,
		resolution = 300
		)

### FIGURE 6G #####################################################################################
# read in data 
setwd('/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_6/')


files <- list.files(path = 'Figure_6_G', pattern = 'markers', full.name = TRUE)
gsea <- sapply(files, function(x) {
	tmp <- read.delim(x, as.is = TRUE)
	tmp[order(rownames(tmp)),]
	}, 
	simplify = FALSE)
plot_data_es <- do.call(cbind, lapply(gsea, '[[', 1))
plot_data_fdr <- do.call(cbind, lapply(gsea, '[[', 3))

plot_data_es

colnames(plot_data_es) <- sapply(colnames(plot_data_es), function(x) paste(strsplit(x, '_')[[1]][4:5], collapse = '_'))
colnames(plot_data_fdr) <- sapply(colnames(plot_data_fdr), function(x) paste(strsplit(x, '_')[[1]][4:5], collapse = '_'))
rownames(plot_data_es) <- rownames(gsea[[1]])
rownames(plot_data_fdr) <- rownames(gsea[[1]])



# only keep a subset of pathways
pathways <- c('HALLMARK_TNFA_SIGNALING_VIA_NFKB','HALLMARK_HYPOXIA','HALLMARK_APOPTOSIS','HALLMARK_CHOLESTEROL_HOMEOSTASIS',
	'HALLMARK_MTORC1_SIGNALING','HALLMARK_P53_PATHWAY','HALLMARK_OXIDATIVE_PHOSPHORYLATION','HALLMARK_MYC_TARGETS_V1',
	'HALLMARK_E2F_TARGETS','HALLMARK_DNA_REPAIR','HALLMARK_G2M_CHECKPOINT')
plot_data_es <- plot_data_es[pathways,]
plot_data_fdr <- plot_data_fdr[pathways,]

# split into two plots
p1samples <- c('D2C2_MID','D2C2_LATE','R1T12_0','R2T12_0','R3T12_0')
plot_data_es_p1 <- plot_data_es[,p1samples]
plot_data_fdr_p1 <- plot_data_fdr[,p1samples]
p2samples <- c('R2T2_0a','R2T2_0','R2T2_1','R2T2_39','R2T2_13', 'R2T2_56',
	'R2T2_11','R2T2_6','R2T2_9')
plot_data_es_p2 <- plot_data_es[,p2samples]
plot_data_fdr_p2 <- plot_data_fdr[,p2samples]

# set yaxis lab 
yaxis.lab <- rownames(plot_data_es)
yaxis.lab <- gsub("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "TNFA", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_HYPOXIA", "HYPOXIA", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_APOPTOSIS", "APOPT.", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_CHOLESTEROL_HOMEOSTASIS", "CHOL HOM", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_MTORC1_SIGNALING", "MTORC1", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_P53_PATHWAY", "P53", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "OX PHOS", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_MYC_TARGETS_V1", "MYC V1", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_E2F_TARGETS", "E2F", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_DNA_REPAIR", "DNA REP", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_G2M_CHECKPOINT", "G2M", yaxis.lab)

# set spot size and colour functions
spot.size.function <- function(x) {abs(x)*3}
spot.colour.function <- function(x) {
    colours <- rep('white', length(x));
    colours[x > 0] <- 'firebrick3';
    colours[x < 0] <- 'dodgerblue3';
    colours[x == 0] <- 'transparent';
    return(colours);
    }

timecov <- c('#CFE1B9','#97A97C','#6C7854','#404A30')
p1cov <- timecov[c(2,4,3,3,3)]
p2cov <- rep(timecov[1],9)

cov1 <- list(
 	rect = list(
 		col = 'transparent',
 		fill = p1cov
 		)
 	);
cov1.grob <- covariates.grob(
	covariates = cov1,
	ord = c(1:length(p1cov)),
	side = 'top',
	size = 1
	);

p1 <- create.dotmap(
		x = plot_data_es_p1,
		bg.data = -log10(plot_data_fdr_p1),
		xaxis.cex = 1.2,
		yaxis.tck = 0,
		xaxis.tck = 0,
		spot.colour.function = spot.colour.function,
		spot.size.function = spot.size.function,
		yaxis.lab = yaxis.lab,
		xaxis.lab = rep('', ncol(plot_data_es_p1)),
		colour.scheme = c('white','black'),
		na.spot.size = 0,
		legend = list(
            bottom = list(fun = cov1.grob)
            ),
		bg.alpha = 1,
		at = c(0, seq(1.3, 6, 0.2))
		)

cov2 <- list(
 	rect = list(
 		col = 'transparent',
 		fill = p2cov
 		)
 	);
cov2.grob <- covariates.grob(
	covariates = cov2,
	ord = c(1:length(p2cov)),
	side = 'top',
	size = 1
	);

p2 <- create.dotmap(
		x = plot_data_es_p2,
		bg.data = -log10(plot_data_fdr_p2),
		xaxis.cex = 1.2,
		yaxis.tck = 0,
		xaxis.tck = 0,
		#xaxis.lab.top = c('0a','0b','1','39','13','56','11','6','9'),
		spot.colour.function = spot.colour.function,
		spot.size.function = spot.size.function,
		yaxis.lab = rep('', nrow(plot_data_es_p2)),
		xaxis.lab = rep('', ncol(plot_data_es_p2)),
		colour.scheme = c('white','black'),
		na.spot.size = 0,
		legend = list(
            bottom = list(fun = cov2.grob)
            ),
		bg.alpha = 1,
		at = c(0, seq(1.3, 6, 0.2))
		)

key_sizes <- seq(1, -1, -0.5);
key_cov <- draw.key(
	list(
			space = 'right',
			points = list(
				cex = spot.size.function(key_sizes),
				col = 'black',
				fill = spot.colour.function(key_sizes),
				pch = 21
				),
			text = list(
				lab = c('1.0','0.5','0.0','-0.5','-1.0'),
				cex = 1.5,
				adj = 1.5,
				fontface = 'bold'
				),
			title = expression(bold("Observed\nScore")),
			cex.title = 1.8,
			background = 'white',
			padding = 8
			)
	)

cov.legend <- list(
	legend = list(
		colours = c('#CFE1B9','#97A97C','#6C7854','#404A30'),
        labels = c('173d','296d','315d','729d'),
		title = expression(bold('Timepoint')),
		cex = 2,
		border = 'transparent'
		)
	);

cov.legend.grob <- legend.grob(
	legends = cov.legend,
	label.cex = 2,
	title.cex = 2
	)

pdf('Figure_6_G/Figure6G_legend.pdf', height = 6.75, width = 12)
grid.arrange(
	key_cov,
	cov.legend.grob,
	ncol = 6
	)
dev.off()

create.multipanelplot(
	list(p1, p2),
	filename = 'Figure_6_GFigure6G.pdf',
	layout.height = 1,
	layout.width = 2,
	legend = list(
			bottom = list(
				fun = BoutrosLab.plotting.general::create.colourkey(
                    x = -log10(plot_data_fdr_p2),
                    colour.scheme = c('white', 'black'),
                    at = c(0, seq(1.3, 6.1, 0.2)),
				    colourkey.labels.at = seq(0,6,2),
					colourkey.labels.cex = 1.5,
                     colourkey.labels = c(
						expression(1),
						expression(10^-2),
						expression(10^-4),
						expression(10^-6)
						),
                    placement = viewport(just = 'left', x = 0.137, y = 0.8, width = 0.86)
                    )
				)
			# inside = list(
			# 		fun = grid.arrange(
			# 			key_cov,
			# 			cov.legend.grob,
			# 			padding = -10,
			# 			ncol = 2
			# 			),
			# 		x = 0.1,
			# 		y = 0.7
			# 		)
			),
	plot.objects.width = c(6.2, 7.8),
	resolution = 300,
	#top.padding = 10,
	height = 6.75,
	width = 12
	)

### EXTENDED DATA 9 #######################################################################
# read in data 
files <- list.files(path = 'Katie_suppFig_S18', pattern = 'markers', full.name = TRUE)
gsea <- sapply(files, function(x) {
	tmp <- read.delim(x, as.is = TRUE)
	tmp[order(rownames(tmp)),]
	}, 
	simplify = FALSE)
plot_data_es <- do.call(cbind, lapply(gsea, '[[', 1))
plot_data_fdr <- do.call(cbind, lapply(gsea, '[[', 3))

colnames(plot_data_es) <- sapply(colnames(plot_data_es), function(x) paste(strsplit(x, '_')[[1]][4:5], collapse = '_'))
colnames(plot_data_fdr) <- sapply(colnames(plot_data_fdr), function(x) paste(strsplit(x, '_')[[1]][4:5], collapse = '_'))
rownames(plot_data_es) <- rownames(gsea[[1]])
rownames(plot_data_fdr) <- rownames(gsea[[1]])

# only keep a subset of pathways
pathways <- c('HALLMARK_TNFA_SIGNALING_VIA_NFKB','HALLMARK_HYPOXIA','HALLMARK_APOPTOSIS','HALLMARK_CHOLESTEROL_HOMEOSTASIS',
	'HALLMARK_MTORC1_SIGNALING','HALLMARK_P53_PATHWAY','HALLMARK_OXIDATIVE_PHOSPHORYLATION','HALLMARK_MYC_TARGETS_V1',
	'HALLMARK_E2F_TARGETS','HALLMARK_DNA_REPAIR','HALLMARK_G2M_CHECKPOINT')
plot_data_es <- plot_data_es[pathways,]
plot_data_fdr <- plot_data_fdr[pathways,]

# split into two plots
p1samples <- c('D3C2_MID','D3C2_LATE')
plot_data_es_p1 <- plot_data_es[,p1samples]
plot_data_fdr_p1 <- plot_data_fdr[,p1samples]
p2samples <- c('R1T2_0a','R1T2_0b','R1T2_2a','R1T2_2b','R1T2_2', 'R1T2_7',
	'R1T2_8','R1T2_3','R1T2_4','R1T2_1')
plot_data_es_p2 <- plot_data_es[,p2samples]
plot_data_fdr_p2 <- plot_data_fdr[,p2samples]

# set yaxis lab 
yaxis.lab <- rownames(plot_data_es)
yaxis.lab <- gsub("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "TNFA", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_HYPOXIA", "HYPOXIA", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_APOPTOSIS", "APOPT.", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_CHOLESTEROL_HOMEOSTASIS", "CHOL HOM", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_MTORC1_SIGNALING", "MTORC1", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_P53_PATHWAY", "P53", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "OX PHOS", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_MYC_TARGETS_V1", "MYC V1", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_E2F_TARGETS", "E2F", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_DNA_REPAIR", "DNA REP", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_G2M_CHECKPOINT", "G2M", yaxis.lab)

# set spot size and colour functions
spot.size.function <- function(x) {abs(x)*3}
spot.colour.function <- function(x) {
    colours <- rep('white', length(x));
    colours[x > 0] <- 'firebrick3';
    colours[x < 0] <- 'dodgerblue3';
    colours[x == 0] <- 'transparent';
    return(colours);
    }


p1 <- create.dotmap(
		x = plot_data_es_p1,
		bg.data = -log10(plot_data_fdr_p1),
		xaxis.cex = 1.2,
		yaxis.tck = 0,
		xaxis.tck = 0,
		spot.colour.function = spot.colour.function,
		spot.size.function = spot.size.function,
		yaxis.lab = yaxis.lab,
		xaxis.lab = rep('', ncol(plot_data_es_p1)),
		colour.scheme = c('white','black'),
		na.spot.size = 0,
		# legend = list(
  #           bottom = list(fun = cov1.grob)
  #           ),
		bg.alpha = 1,
		at = c(0, seq(1.3, 6.1, 0.2))
		)


p2 <- create.dotmap(
		x = plot_data_es_p2,
		bg.data = -log10(plot_data_fdr_p2),
		xaxis.cex = 1.2,
		yaxis.tck = 0,
		xaxis.tck = 0,
		#xaxis.lab.top = c('0a','0b','1','39','13','56','11','6','9'),
		spot.colour.function = spot.colour.function,
		spot.size.function = spot.size.function,
		yaxis.lab = rep('', nrow(plot_data_es_p2)),
		xaxis.lab = rep('', ncol(plot_data_es_p2)),
		colour.scheme = c('white','black'),
		na.spot.size = 0,
		# legend = list(
  #           bottom = list(fun = cov2.grob)
  #           ),
		bg.alpha = 1,
		at = c(0, seq(1.3, 6.1, 0.2))
		)

key_sizes <- seq(1, -1, -0.5);
key_cov <- draw.key(
	list(
			space = 'right',
			points = list(
				cex = spot.size.function(key_sizes),
				col = 'black',
				fill = spot.colour.function(key_sizes),
				pch = 21
				),
			text = list(
				lab = c('1.0','0.5','0.0','-0.5','-1.0'),
				cex = 1.5,
				adj = 1.5,
				fontface = 'bold'
				),
			title = expression(bold("Observed\nScore")),
			cex.title = 1.8,
			background = 'white',
			padding = 8
			)
	)

pdf('FigureS18_legend.pdf', height = 6.75, width = 12)
grid.draw(key_cov)
dev.off()

create.multipanelplot(
	list(p1, p2),
	filename = 'FigureS18.pdf',
	layout.height = 1,
	layout.width = 2,
	legend = list(
			bottom = list(
				fun = BoutrosLab.plotting.general::create.colourkey(
                    x = -log10(plot_data_fdr_p2),
                    colour.scheme = c('white', 'black'),
                    at = c(0, seq(1.3, 6.1, 0.2)),
				    colourkey.labels.at = seq(0,6,2),
					colourkey.labels.cex = 1.5,
                     colourkey.labels = c(
						expression(1),
						expression(10^-2),
						expression(10^-4),
						expression(10^-6)
						),
                    placement = viewport(just = 'left', x = 0.137, y = 0.8, width = 0.86)
                    )
				)
			# inside = list(
			# 		fun = grid.arrange(
			# 			key_cov,
			# 			cov.legend.grob,
			# 			padding = -10,
			# 			ncol = 2
			# 			),
			# 		x = 0.1,
			# 		y = 0.7
			# 		)
			),
	plot.objects.width = c(2.8, 7.2),
	resolution = 300,
	#top.padding = 10,
	height = 6.75,
	width = 12
	)

### EXTENDED DATA 7 #######################################################################
# read in data 
files <- list.files(path = 'Katie_suppFig_S16', pattern = 'markers', full.name = TRUE)
gsea <- sapply(files, function(x) {
	tmp <- read.delim(x, as.is = TRUE)
	tmp[order(rownames(tmp)),]
	}, 
	simplify = FALSE)
plot_data_es <- do.call(cbind, lapply(gsea, '[[', 1))
plot_data_fdr <- do.call(cbind, lapply(gsea, '[[', 3))

colnames(plot_data_es) <- sapply(colnames(plot_data_es), function(x) paste(strsplit(x, '_')[[1]][4:5], collapse = '_'))
colnames(plot_data_fdr) <- sapply(colnames(plot_data_fdr), function(x) paste(strsplit(x, '_')[[1]][4:5], collapse = '_'))
rownames(plot_data_es) <- rownames(gsea[[1]])
rownames(plot_data_fdr) <- rownames(gsea[[1]])

# only keep a subset of pathways
pathways <- c('HALLMARK_TNFA_SIGNALING_VIA_NFKB','HALLMARK_HYPOXIA','HALLMARK_APOPTOSIS','HALLMARK_CHOLESTEROL_HOMEOSTASIS',
	'HALLMARK_MTORC1_SIGNALING','HALLMARK_P53_PATHWAY','HALLMARK_OXIDATIVE_PHOSPHORYLATION','HALLMARK_MYC_TARGETS_V1',
	'HALLMARK_E2F_TARGETS','HALLMARK_DNA_REPAIR','HALLMARK_G2M_CHECKPOINT')
plot_data_es <- plot_data_es[pathways,]
plot_data_fdr <- plot_data_fdr[pathways,]

# set spot size and colour functions
spot.size.function <- function(x) {abs(x)*3}
spot.colour.function <- function(x) {
    colours <- rep('white', length(x));
    colours[x > 0] <- 'firebrick3';
    colours[x < 0] <- 'dodgerblue3';
    colours[x == 0] <- 'transparent';
    return(colours);
    }

key_sizes <- seq(1, -1, -0.5);
key_cov <- draw.key(
	list(
			space = 'right',
			points = list(
				cex = spot.size.function(key_sizes),
				col = 'black',
				fill = spot.colour.function(key_sizes),
				pch = 21
				),
			text = list(
				lab = c('1.0','0.5','0.0','-0.5','-1.0'),
				cex = 1.5,
				adj = 1.5,
				fontface = 'bold'
				),
			title = expression(bold("Observed\nScore")),
			cex.title = 1.8,
			background = 'white',
			padding = 8
			)
	)

create.dotmap(
		x = plot_data_es,
		bg.data = -log10(plot_data_fdr),
		filename = 'FigureS16.pdf',
		xaxis.cex = 1.2,
		yaxis.tck = 0,
		xaxis.tck = 0,
		spot.colour.function = spot.colour.function,
		spot.size.function = spot.size.function,
		yaxis.lab = yaxis.lab,
		xaxis.lab = rep('', ncol(plot_data_es)),
		colour.scheme = c('white','black'),
		colourkey = TRUE,
		colourkey.labels.at = seq(0,6,2),
		colourkey.cex = 1.5,
		colourkey.labels = c(
						expression(1),
						expression(10^-2),
						expression(10^-4),
						expression(10^-6)
						),
		key = list(
			space = 'left',
			points = list(
				cex = spot.size.function(key_sizes),
				col = 'black',
				fill = spot.colour.function(key_sizes),
				pch = 21
				),
			text = list(
				lab = c('1.0','0.5','0.0','-0.5','-1.0'),
				cex = 1.5,
				adj = 1.5,
				fontface = 'bold'
				),
			title = expression(bold("Observed\nScore")),
			cex.title = 1.8,
			background = 'white',
			padding = 8
			),
		resolution = 300,
		height = 6.75,
		width = 10.5,
		right.padding = 2,
		bottom.padding = 5,
		left.padding = 3,
		na.spot.size = 0,
		bg.alpha = 1,
		at = c(0, seq(1.3, 6.1, 0.2))
		)

### EXTENDED DATA 8 #######################################################################
# read in data 
files <- list.files(path = 'Katie_suppFig_S17', pattern = 'markers', full.name = TRUE)
gsea <- sapply(files, function(x) {
	tmp <- read.delim(x, as.is = TRUE)
	tmp[order(rownames(tmp)),]
	}, 
	simplify = FALSE)
plot_data_es <- do.call(cbind, lapply(gsea, '[[', 1))
plot_data_fdr <- do.call(cbind, lapply(gsea, '[[', 3))

colnames(plot_data_es) <- sapply(colnames(plot_data_es), function(x) paste(strsplit(x, '_')[[1]][4:5], collapse = '_'))
colnames(plot_data_fdr) <- sapply(colnames(plot_data_fdr), function(x) paste(strsplit(x, '_')[[1]][4:5], collapse = '_'))
rownames(plot_data_es) <- rownames(gsea[[1]])
rownames(plot_data_fdr) <- rownames(gsea[[1]])

# only keep a subset of pathways
pathways <- c('HALLMARK_TNFA_SIGNALING_VIA_NFKB','HALLMARK_HYPOXIA','HALLMARK_APOPTOSIS','HALLMARK_CHOLESTEROL_HOMEOSTASIS',
	'HALLMARK_MTORC1_SIGNALING','HALLMARK_P53_PATHWAY','HALLMARK_OXIDATIVE_PHOSPHORYLATION','HALLMARK_MYC_TARGETS_V1',
	'HALLMARK_E2F_TARGETS','HALLMARK_DNA_REPAIR','HALLMARK_G2M_CHECKPOINT')
plot_data_es <- plot_data_es[pathways,]
plot_data_fdr <- plot_data_fdr[pathways,]

# set spot size and colour functions
spot.size.function <- function(x) {abs(x)*3}
spot.colour.function <- function(x) {
    colours <- rep('white', length(x));
    colours[x > 0] <- 'firebrick3';
    colours[x < 0] <- 'dodgerblue3';
    colours[x == 0] <- 'transparent';
    return(colours);
    }

key_sizes <- seq(1, -1, -0.5);
key_cov <- draw.key(
	list(
			space = 'right',
			points = list(
				cex = spot.size.function(key_sizes),
				col = 'black',
				fill = spot.colour.function(key_sizes),
				pch = 21
				),
			text = list(
				lab = c('1.0','0.5','0.0','-0.5','-1.0'),
				cex = 1.5,
				adj = 1.5,
				fontface = 'bold'
				),
			title = expression(bold("Observed\nScore")),
			cex.title = 1.8,
			background = 'white',
			padding = 8
			)
	)

create.dotmap(
		x = plot_data_es,
		bg.data = -log10(plot_data_fdr),
		filename = 'FigureS17.pdf',
		xaxis.cex = 1.2,
		yaxis.tck = 0,
		xaxis.tck = 0,
		spot.colour.function = spot.colour.function,
		spot.size.function = spot.size.function,
		yaxis.lab = yaxis.lab,
		xaxis.lab = rep('', ncol(plot_data_es)),
		colour.scheme = c('white','black'),
		colourkey = TRUE,
		colourkey.labels.at = seq(0,6,2),
		colourkey.cex = 1.5,
		colourkey.labels = c(
						expression(1),
						expression(10^-2),
						expression(10^-4),
						expression(10^-6)
						),
		key = list(
			space = 'left',
			points = list(
				cex = spot.size.function(key_sizes),
				col = 'black',
				fill = spot.colour.function(key_sizes),
				pch = 21
				),
			text = list(
				lab = c('1.0','0.5','0.0','-0.5','-1.0'),
				cex = 1.5,
				adj = 1.5,
				fontface = 'bold'
				),
			title = expression(bold("Observed\nScore")),
			cex.title = 1.8,
			background = 'white',
			padding = 8
			),
		resolution = 300,
		height = 6.75,
		width = 8.5,
		right.padding = 2,
		bottom.padding = 5,
		left.padding = 3,
		na.spot.size = 0,
		bg.alpha = 1,
		at = c(0, seq(1.3, 6.1, 0.2))
		)

#### SEXTENDED DATA 10 #####################################################################
# read in gsea results
pathwayorder <- read.delim(
	"Altered_pathway_counts_order.txt",
	header=TRUE
	)
all <- read.csv(
	"GSEA_compiled_all_direction_full_combined2.csv",
	header=TRUE
	)
all$p_val_adj <- abs(all$p_val_adj)
all[all$p_val_adj > 0.05,'Observed_Score'] <- 0
# reformat data
plot_data_es <- spread(
	all[,c('Pathway','Sample','Observed_Score')],
	key = Sample,
	value = Observed_Score
	)
plot_data_fdr <- spread(
	all[,c('Pathway','Sample','p_val_adj')],
	key = Sample,
	value = p_val_adj
	)

### Order pathways based on the exact fisher results
plot_data_es$Pathway <- factor(plot_data_es$Pathway, levels = pathwayorder$Pathway)
plot_data_fdr$Pathway <- factor(plot_data_fdr$Pathway, levels = pathwayorder$Pathway)
plot_data_es <- plot_data_es[order(plot_data_es$Pathway),]
plot_data_fdr <- plot_data_fdr[order(plot_data_fdr$Pathway),]

plot_data_es$Pathway <- gsub('HALLMARK_', '', plot_data_es$Pathway)
plot_data_fdr$Pathway <- gsub('HALLMARK_', '', plot_data_fdr$Pathway)
rownames(plot_data_es) <- plot_data_es$Pathway
rownames(plot_data_fdr) <- plot_data_fdr$Pathway
plot_data_es <- plot_data_es[,-1]
plot_data_fdr <- plot_data_fdr[,-1]

### Order samples
samples <- c("D1C1_LATE","D1C2_LATE","D1C3_LATE","D2C2_LATE","D2C3_LATE","D3C2_LATE",
                       "D2C1R1_0a","D2C1R2_1a","D2C2R2_0a","D3C2R1_0a",
                       "D2C1R1_0b","D2C1R1_1","D2C1R1_2","D2C1R1_4",
                       "D2C1R2_0a","D2C1R2_0b","D2C1R2_1b","D2C1R2_2","D2C1R2_4",
                       "D2C2R2_0b","D2C2R2_1","D2C2R2_11","D2C2R2_13","D2C2R2_39","D2C2R2_56","D2C2R2_6","D2C2R2_9",
                       "D3C2R1_0b","D3C2R1_1","D3C2R1_2a","D3C2R1_2b","D3C2R1_2c","D3C2R1_3","D3C2R1_4","D3C2R1_7","D3C2R1_8")
plot_data_es <- plot_data_es[,samples]
plot_data_fdr <- plot_data_fdr[,samples]

# set yaxis lab 
yaxis.lab <- rownames(plot_data_es)
yaxis.lab <- gsub("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "TNFA", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_HYPOXIA", "HYPOXIA", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_APOPTOSIS", "APOPT.", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_CHOLESTEROL_HOMEOSTASIS", "CHOL HOM", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_MTORC1_SIGNALING", "MTORC1", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_P53_PATHWAY", "P53", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "OX PHOS", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_MYC_TARGETS_V1", "MYC V1", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_E2F_TARGETS", "E2F", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_DNA_REPAIR", "DNA REP", yaxis.lab)
yaxis.lab <- gsub("HALLMARK_G2M_CHECKPOINT", "G2M", yaxis.lab)

# set spot size and colour functions
spot.size.function <- function(x) {abs(x)*2}
spot.colour.function <- function(x) {
    colours <- rep('white', length(x));
    colours[x > 0] <- 'firebrick3';
    colours[x < 0] <- 'dodgerblue3';
    colours[x == 0] <- 'transparent';
    return(colours);
    }

# set covariate to show clusters and cohort
sample_cov <- c(rep('#772766', 6), rep('#FF0000', 4), rep('#D1D3D4', 26))

top_covariate <- list(
	rect = list(
		col = 'white',
		fill = sample_cov,
		lwd = 1.5
		)
	)
top_cov_grob <- covariates.grob(
	covariates = top_covariate,
	ord = 1:length(sample_cov),
	side = 'top'
	)

right_covariate <- list(
	rect = list(
		col = 'white',
		fill = c(rep('#F389AA', 10), rep('white', 40)),
		lwd = 1.5
		)
	)
right_cov_grob <- covariates.grob(
	covariates = right_covariate,
	ord = 50:1,
	side = 'right'
	)


cov_legend <- list(
	legend = list(
		colours = "#F389AA",
		labels = 'Top Altered\nPathways',
		title = ''
		),
	legend = list(
		colours = c("#772766", "#FF0000","#D1D3D4"),
		labels = c('Late\nculture','Winning\nsubclones','Losing\nsubclones'),
		title = ''
		)
	)

legend_grob <- legend.grob(
	legends = cov_legend,
	label.cex = 1.8,
	title.cex = 1.8
	)

key_sizes <- seq(1, -1, -0.5);
key_cov <- draw.key(
	list(
			space = 'right',
			points = list(
				cex = spot.size.function(key_sizes),
				col = 'black',
				fill = spot.colour.function(key_sizes),
				pch = 21
				),
			text = list(
				lab = c('1.0','0.5','0.0','-0.5','-1.0'),
				cex = 1.5,
				adj = 1.5,
				fontface = 'bold'
				),
			title = expression(bold("Effect\nSize")),
			cex.title = 1.8,
			background = 'white',
			padding = 8
			)
	)

# create dotmap
create.dotmap(
		x = plot_data_es,
		bg.data = -log10(plot_data_fdr),
		filename = 'Figure22A.pdf',
		xaxis.cex = 1,
		yaxis.cex = 1,
		yaxis.tck = 0,
		xaxis.tck = 0,
		spot.colour.function = spot.colour.function,
		spot.size.function = spot.size.function,
		yaxis.lab = yaxis.lab,
		xaxis.lab = paste('     ', colnames(plot_data_es)),
		xaxis.rot = 90,
		colour.scheme = c('white','black'),
		na.spot.size = 0,
		row.colour = 'white',
		col.colour = 'white',
		bg.alpha = 1,
		colourkey = TRUE,
		colourkey.labels = c(
			expression(1),
			expression(10^-2),
			expression(10^-4),
			expression(10^-6)
			),
		at = c(0, seq(1.3, 6.1, 0.2)),
		colourkey.labels.at = seq(0,6,2),
		colourkey.cex = 1.5,
		axis.top = 1.5,
		key = NULL,
		legend = list(
			right = list(
				fun = right_cov_grob
				),
			top = list(
				fun = top_cov_grob
				),
			inside = list(
				fun = legend_grob,
				x = 1.1,
				y = 1.03
				),
			inside = list(
				fun = key_cov,
				x = 1.1,
				y = 0.55
				),
			inside = list(
				fun = draw.key(list(text = list(lab = ''), title = expression(bold('FDR')), cex.title = 1.8)),
				x = -0.1,
				y = -0.125
				)
			),
		width = 13,
		height = 12,
		bottom.padding = 5,
		right.padding = 27,
		top.padding = 5,
		left.padding = 5,
		resolution = 300
		)


