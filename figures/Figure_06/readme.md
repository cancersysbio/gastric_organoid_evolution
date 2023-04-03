# ECB cDNA analysis

**Figure_6A.1_inferCNV_prepare_Seq8_input.R** - R script that takes the count matrix from Seurat, after the removal of bad quality cells and wrangles it into shape for the analysis script. Here, data from D2C2 replicate 2 at day 168 is used specifically.

**Figure_6A.2_inferCNV_analysis.R** - R script with the commands to run inferCNV.

**Figure_6A.3_sbatch_Seq8_inferCNV.sh** - Bash script to run inferCNV on the data from D2R2 replicate 2 at day 173 on a slurm cluster.

**Figure_6A.4_Seq8_ECB_heatmaps.R** - Script to re-plot the inferCNV output with cells ordered according to their cell-subclone-relationship. The resulting heatmap is shown in Figure 6a.

**Figure_6A.5_Seq8_ECB_heatmap_chr.R** - Script to plot a zoom-in of specific chromosomes from D2C2 replicate 2 day 168, that show copy number alterations between the winning subclone and the base-line subclone. The zoom-in is represented in Figure 6a.

**Figure_6C.1_fishplot_barcode_frequencies.py** - Rescaling and reformatting for fishplot. Used for figure 6b

**Figure_6C.2_barcode_fishplot.R** - Script for combined CNA and barcode fishplot. Used for figure 6b

**Figure_6D_growth_vs_frequency.R** - Script to plot barcode growth versus frequency figure. Used for figure 6c

**Figure_6E_expression_dotplot.R** - Seurat analysis of ECB data. Used for figures 6d-g

**Figure_6F_volcano_plot.R** - Script for volcano plots. Used for figure 6e

**Figure_6H.1_Compile_list_of_GSEA_direction.py** - Script to compile ECB GSEA output

**Figure_6H.2_Compile_list_of_GSEA_EML_direction.py** - Script to compile EML GSEA output

**Figure_6H.3_GSEA_correlations.R** - Script to calculate and plot GSEA sample to sample correlations. Used for figure 6g

**Note:**

**Figure_6G** - Scripts for GSEA analysis and plotting are the same as the ones found under figure 3
