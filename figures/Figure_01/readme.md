# Shallow whole genome sequencing analysis

## Scripts

**Figure_1B-C_sWGS_CNA_plot.R** - Script to generate scatterplots from shallow wgs results processed with qdnaseq. Source data 50k sWGS *CNA_smooth.txt, output from qDNAseq. Source data files found on Zenodo

**Figure_1D.1_chr_smooth_add_chr_arm.py** - Script to amend qdnaseq output files "cnv_smoot" with chromosome arm location, used for both figure 1 and 2 (sWGS and WGS)

**Figure_1D.2_moving_average.py** - Script to add moving average to qdnaseq output files and transforms reads per window to read counts centered around 1 (instead of log 2 value), used for both figure 1 and 2 (shallow and deep sequencing)

**Figure_1D.3_fraction_genome_altered.py** - Script to summarize fraction genome altered per culture, used in figure 1d

**Figure_1D.4_plot_fraction_genome_altered.R** - Script to generate scatterplot for fraction genome altered

**Figure_1E.1_calculate_biallelic_loss.py** - Script to calculate segments of bi-allelic loss, egments are called as lost if a genomic window contains 25% or less reads compared to median number of reads per window

**Figure_1E.2_biallelic_loss_first_occurrance.py** - Script to calculate biallelic loss, the program outputs the first occurance of a biallelic loss

**Figure_1E.3_aberration_timing.py** - Script to calculate genome aberration timing, outputs for each sample if a chromosome arm has a copy number aberration or not

**Figure_1E.4_aberration_timing_first_occurrance.py** - Script to calculate the timing of copy number alterations, the program outputs the first occurrance of a copy number aberration

**Figure_1E.5_plot_timing_genomic_events.R** - Script to generate boxplot for timing of genomic events
