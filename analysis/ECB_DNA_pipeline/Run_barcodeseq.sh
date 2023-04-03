#!/bin/bash

### set the name of the job
#SBATCH --job-name=runBarcodeSeq

### Job output files
#SBATCH --output=ecb.%j.out
#SBATCH --error=ecb.%j.err

### Set working directory
#SBATCH --chdir=/N/users/kasperk/storage/ECB/Longitudinal_analysis/sequtils/ECB_newcode

### Set job time
#SBATCH --time=0-08:00:00

### Allocate resources
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G

### Set notifications

#SBATCH --mail-type=END,FAIL # notifications for job done & fail

### Start your program
### Python Settings

### SET DIRECTORY CONTAINING FASTQ FILES AND THE SUBDIRECTORY
#dir="/N/users/kasperk/storage/ECB/Longitudinal_analysis/samples/BarcodeSeq/D1C1"
#dir="/N/users/kasperk/storage/ECB/Longitudinal_analysis/samples/BarcodeSeq/D2C1"
#dir="/N/users/kasperk/storage/ECB/Longitudinal_analysis/samples/BarcodeSeq/D2C2"
dir="/N/users/kasperk/storage/ECB/Longitudinal_analysis/samples/BarcodeSeq/D2C3"
#dir="/N/users/kasperk/storage/ECB/Longitudinal_analysis/samples/BarcodeSeq/D3C2"

### SET SUBDIRECTORY CONTAINING NAME FILE OF SAMPLE TO USE
#subdir="T8"
#subdir="T12"
subdir="T20"

### SPECIFY NAME OF FILE WITH SAMPLES YOU WANT TO USE
#namefile=samples_D1C1_T12
#namefile=samples_D2C1_T20
#namefile=samples_D2C1_T12_correct
#namefile=samples_T12_D2C2_correct
namefile=samples_T20_D2C2_correct
#namefile=samples_T12_D2C3_correct
namefile=samples_T20_D2C3_correct
#namefile=samples_D3C2_T20_correct
#namefile=samples_D3C2_T12_correct

### SET NAME OF ECB MASTER LIST AND RG OUTPUT FOLDER
### BASE
#masterfile="ECB_groups_master_list.txt"  ### D1C1, D2C3_T12, D3C2_T12, D2C1_T12
masterfile="ECB_groups_master_list_T12.txt"  ### D2C3_T20
RG_folder="RG_base"



### PAIRS
#masterfile="ECB_groups_master_list_updated.txt"  ### D2C1, D2C2, D3C2
#RG_folder="RG_pairs"


### R Settings
Rdir=$dir/$subdir/$RG_folder
#nrParents="1" # For D1C1 and D3C2
nrParents="2" # For D2C1, D2C2 and D2C3


pythonpath="/N/users/kasperk/anaconda3/bin/python"

$pythonpath barcodeseq_parsing.py $dir $subdir $namefile $masterfile $RG_folder

source ~/.bashrc
conda activate r_env_seurat

Rscript BarcodeSeq_plotting.R $Rdir $nrParents

