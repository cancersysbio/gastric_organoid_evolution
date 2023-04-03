# ECB cDNA pipeline

**ecb_functions.py** - Functions to extract ECB cDNA barcodes reads from fastq files

**ecb_parsing.py** - Functions to associate ECB cDNA barcodes with prefetermined ECB read groups and associate these with specific cells in the single cell RNA sequencing data. Also functions to call if multiple ECB tags were inserted during transduction.

**hash_functions.py** - Functions to extract cell hashing barcode from fastq files

**hash_parsing.py** - Functions to associate hash barcodes with specific cells in the single cell RNA sequencing data

**combine_ecb_hash_bcrg_updated.py** - Functions to combine the ECB and hash outputs. "updated" refers to that calling of multiple ECB insertions during transductions should be done in this sample and that output files should be updated

**combine_ecb_hash_bcrg_not_updated.py** - Functions to combine the ECB and hash outputs. "not updated" refers to that calling of multiple ECB insertions during transductions has been done on another samples, and that output files should not be updated based on this calling.

**Bash_forECBandHash.sh** - Shell script to run ecb_parsing.py, hash_parsing.py, combine_ecb_hash_bcrg_updated.py and combine_ecb_hash_bcrg_not_updated.py
