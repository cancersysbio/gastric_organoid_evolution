# Quantifying evolution towards malignancy in Early and Late organoids

**Figure_4A-B.1_Seurat_SatheEtal.R** - Script to re-analyze the publicly available dataset from Sathe et al., 2020, Clinical Cancer Research using Seurat (<https://clincancerres.aacrjournals.org/content/clincanres/early/2020/04/03/1078-0432.CCR-19-3231.full.pdf>).

**Figure_4A-B.2_COMBAT.R** - Script for the removal of the organoid-specific gene signature by using reference data from Sathe et al.

**Figure_4A-B.3_Clustering_UMAP_Sathe_EML.R** - Script for clustering of the scRNA-seq data from Sathe et al., generated with *Figure_4A-B.1_Seurat_SatheEtal.R*, which are then subsequently visualized using UMAP. The UMAP representation is then saved for the projection in *Figure_4C-E_Sathe_EML.R*. This script generates Figure 4a and b.

**Figure_4C-E_Sathe_EML.R** - Here, the LSI projection itself is implemented. This script projects each individual organoid clone at the *Early*, *Mid* and *Late* time point into the gastric dataset form *Figure_4A-B.2_COMBAT.R*. Subsequently, the nearest neighbors are quantified and used to quantify the NN cell type frequency. This script generates Figure 4c and e.

The scripts represented here were adapted from Jeffrey Granja's [MPAL publication](https://github.com/GreenleafLab/MPAL-Single-Cell-2019).
