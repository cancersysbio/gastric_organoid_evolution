# Whole genome sequencing analysis

## Scripts

**Figure_2A.1_combine_mutations.py** - Script to combine mutations from all .maf files for figure 2a

**Figure_2A.2_WGS_extract_information.py** - Script to get nr SNV, Coding SVN, SNV in Cosmic, Indels and Location for figure 2a

**Figure_2A.3_Kataegis_analysis.R** - 

**Figure_2A.4_wgii_and_lohfrac.R** - Script to get weighted genome instability index  and fraction loss of heterozygosity

**Figure_2A.5_create_input_sigs.R** - Script to create input objects for hdp from mutations.

**Figure_2A.6_run_hdp_multichain.sh** - Script to create run hdp process for component identification.

**Figure_2A.7_create_multiChain_hdp_object.R** - Script to perform signature extraction based on COSMIC and extracted hdp components.

**Figure_2A.8_plot.R** -  Script to plot figure 2A barplots

**Figure_2B_circos_plots.R** - Script to plot the circos plots for figure 2b

**Figure_2C-D_WGS_plots.R** - Script to generate scatterplots from wgs results processed with qdnaseq. Used for figures 2d-e

**Figure_2C_Jabba_plot.R** - Script to plot Jabba figures for figure 2c

**Figure_2E.1_extract_CNV_clones.py** - Script to extract copy number alteration subclones. the program inputs median centered and moving averaged smoothened qdnaseq output files and outputs regions that are at least 100 windows long and where each window has at least a 12.5% read count change compared to median read count. Used for figure 2e

**Figure_2E.2_fishplot.R** - Script to generate fishplot of CNA-based subclones in figure 2e
