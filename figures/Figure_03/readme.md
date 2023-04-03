# scRNA-seq

## Data

## Scripts

**Figure_3B.1_inferred_growth_trajectories.R** - Script to calculate cumulative population size over time from cell counts and fit curve to estimate growth derivatives. Used for figure 3b

**Figure_3B.2_inferred_growth_trajectories_make_long.py** - Script to combine growth derivatives, passaging numbers and growth fold change in long format for plotting. Used for figure 3b

**Figure_3B.3_inferred_growth_trajectories_plotting.R** - Script to plot inferred growth derivative in figure 3b

**Figure_3C_3F_scRNA_Seurat_D1.R** - Script to analyze and process the single cell RNA sequencing data from Donor 1 generated from *Early*, *Mid* and *Late* organoids. This script performs the whole Seurat analysis workflow, differential gene expression analysis as well as the plotting panels c and f.

**Figure_3E.1_Seurat_differential_gene_exp_D1.R** - Script to generate differential gene expression for donor 1

**Figure_3E.2_Seurat_differential_gene_exp_D1.R** - Script to generate differential gene expression for donor 2

**Figure_3E.3_Seurat_differential_gene_exp_D1.R** - Script to generate differential gene expression for donor 3

**Figure_3E.4_upset_plot.Rmd** - Script to generate and visualise the upset plots, used for figure 3e

**Figure_3F.1_GSEA_analysis.R** - Script to analyze gene set enrichment from differential expression, used for figures 3j and 6f

**Figure_3F.2_GSEA_plotting.R** - Script to plot gene set enrichment data, used for figures 3j and 6f
