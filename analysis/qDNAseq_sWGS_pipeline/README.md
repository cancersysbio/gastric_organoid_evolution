### QDNAseq copy nummber analysis for sWGS data

We used [`QDNAseq v1.22.0`](https://bioconductor.org/packages/release/bioc/html/QDNAseq.html) for to call copy number aberrations using the low-pass whole genome sequencing data.


### QDNAseqbin annotations

We generated QDNAseq bin annotation for the human genome build hg38 and shared it as a separate [R package](https://github.com/asntech/QDNAseq.hg38)

This package provides QDNAseq bin annotations of size `1, 5, 10, 15, 30, 50, 100, 500 and 1000` kbp for the human genome build hg38.The bin annotations are created using the steps mentioned in QDNAseq vignette and also [here](https://github.com/ccagc/QDNAseq/issues/59).


## Install bin annotations

You need to install the bin annotations v1.1.0 from GitHub:

``` r
#Install the QDNAseq.hg38 package using remotes
remotes::install_github("asntech/QDNAseq.hg38@main")
#or devtools
devtools::install_github("asntech/QDNAseq.hg38@main")
```

### how to run QDNAseq script?

**Please make sure to install the following versions of R/Bioconductor packages before running `qdnaseq.R`.**

```
optparse v1.6.6
tidyverse v1.3.0
Biobase v2.46.0
QDNAseq v1.22.0
CGHbase v1.46.0
future v1.20.1
data.table v1.14.2
ACE v1.4.0
png v0.1.7
grid v3.6.0
ggplot2 v3.3.2
gridExtra v2.3

```

**To run the qdnaseq pipeline use this**

``` shell 
qdnaseq.R --file sample1.bam --out ~/output/dir/ --name sample_name --genome hg38 --minmapq 37 --bin 50
```

To find more details about the script options type `--help`. 
```shell 
qdnaseq.R --help

Options:
	-f CHARACTER, --file=CHARACTER
		Input BAM file or path

	-o CHARACTER, --out=CHARACTER
		output path name [default= ./]

	-n CHARACTER, --name=CHARACTER
		sample bam file name/id [default= sample1]

	-g CHARACTER, --genome=CHARACTER
		Genome (hg19, hg38) [default= hg38]

	-m INTEGER, --minmapq=INTEGER
		Minimum quality score [default= 37]

	-r LOGICAL, --refit=LOGICAL
		Refit the calls using ACE package [default= TRUE]

	-d LOGICAL, --duplicated=LOGICAL
		Are the reads dublicated? [default= FALSE]

	-b INTEGER, --bin=INTEGER
		bin size (kb). Use 5, 10, 15, 30, 50, 100, 200, 500, 1000 [default= 50]

	-h, --help
		Show this help message and exit
```
