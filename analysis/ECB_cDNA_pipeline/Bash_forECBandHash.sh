#!/bin/bash
# To Run: sbatch process_10X_ecb.sh

### Replace MY_PI_SUNetID_or_Project_ID with the PI/project to be charged.
#SBATCH --account=ccurtis2

### set the name of the job
#SBATCH --job-name=ECB_HASH

### Job output files
#SBATCH --output=AAA_ECB_HASH.%j.out
#SBATCH --error=AAA_ECB_HASH.%j.err

### Set working directory
#SBATCH --chdir=/N/users/kasperk/storage/10X/hashECB_data_updated_donors/sequtils

### Set job time
#SBATCH --time=0-08:00:00

### Allocate resources
#SBATCH --cpus-per-task=1
#SBATCH --mem=25G

### Set notifications
#SBATCH --mail-user=kasperk@stanford.edu
#SBATCH --mail-type=END,FAIL # notifications for job done & fail

### Start your program

## Cd to programs paths ##
# make sure it contains: ecb_parsing.py, ecb_functions.py, hash_parsing.py, hash_functions.py
program_path=/N/users/kasperk/storage/10X/hashECB_data_updated_donors/sequtils/
cd ${program_path}


## paths and fastq ##
# ecb_fqpath contains the fq, ecb ref, bc ref, and rg color ref
#ecb_fqpath=/N/users/kasperk/storage/10X/hashECB_data_updated_donors/Sequencing8_D2C2_R2_T2/ECB_data/
ecb_fqpath=/N/users/kasperk/storage/10X/hashECB_data_updated_donors/Sequencing18_D3C2_ECB_T2/ECB_data/
#ecb_fqpath=/N/users/kasperk/storage/10X/hashECB_data_updated_donors/Sequencing20_D2C1_ECB_T2_T7_T12/ECB_data/
#ecb_fqpath=/N/users/kasperk/storage/10X/hashECB_data_updated_donors/Sequencing21_D2C2_ECB_T6/ECB_data/
#ecb_fqpath=/N/users/kasperk/storage/10X/hashECB_data_updated_donors/Sequencing24_D2C1_ECB_R3_T0_T12/ECB_data/
#ecb_fqpath=/N/users/kasperk/storage/10X/hashECB_data_updated_donors/Sequencing25_D2_ECB_T2_T12/ECB_data/

# hash_fqpath contains the fq, adt ref and bc ref
#hash_fqpath=/N/users/kasperk/storage/10X/hashECB_data_updated_donors/Sequencing5_D3C2_thaw/Hash/
#hash_fqpath=/N/users/kasperk/storage/10X/hashECB_data_updated_donors/Sequencing8_D2C2_R2_T2/Hash/
#hash_fqpath=/N/users/kasperk/storage/10X/hashECB_data_updated_donors/Sequencing13_D2_ML/Hash/
hash_fqpath=/N/users/kasperk/storage/10X/hashECB_data_updated_donors/Sequencing18_D3C2_ECB_T2/Hash/
#hash_fqpath=/N/users/kasperk/storage/10X/hashECB_data_updated_donors/Sequencing19_D3_ML/Hash/
#hash_fqpath=/N/users/kasperk/storage/10X/hashECB_data_updated_donors/Sequencing20_D2C1_ECB_T2_T7_T12/Hash/
#hash_fqpath=/N/users/kasperk/storage/10X/hashECB_data_updated_donors/Sequencing21_D2C2_ECB_T6/Hash/
#hash_fqpath=/N/users/kasperk/storage/10X/hashECB_data_updated_donors/Sequencing24_D2C1_ECB_R3_T0_T12/Hash/
#hash_fqpath=/N/users/kasperk/storage/10X/hashECB_data_updated_donors/Sequencing25_D2_ECB_T2_T12/Hash/
#hash_fqpath=/N/users/kasperk/storage/10X/hashECB_data_updated_donors/Sequencing27_D1_ML/Hash/
#hash_fqpath=/N/users/kasperk/storage/10X/hashECB_data_updated_donors/Sequencing29_D1_D2_EARLY/Hash/
#hash_fqpath=/N/users/kasperk/storage/10X/hashECB_data_updated_donors/Sequencing30_D1_D3_EARLY/Hash/
#hash_fqpath=/N/users/kasperk/storage/10X/hashECB_data_updated_donors/Sequencing31_D3C3_EARLY/Hash/

# Location of the input files for combine_ecb_hash.py
ecb_out=${ecb_fqpath}ecb_outs/ # contains BC_RG_AssociatedUMICount_updated.txt
hash_out=${hash_fqpath}hash_outs/intermediate_output/ # contains Hash_correlation.txt

# fastq format: .txt.gz or .fq.gz
# single sample #

#ecb_fq=Seq8_ECB_MOCK_
ecb_fq=Seq18_ECB_TCCGGAGA_
#ecb_fq=Seq20_ECB_ATTACTCG_,Seq20_ECB_TCCGGAGA_
#ecb_fq=Seq21_ECB_ATTACTCG_
#ecb_fq=Seq24_ECB_
#ecb_fq=Seq25_ECB_


#hash_fq=Seq5_Hash_ATTACTCG_,Seq5_Hash_TCCGGAGA_
#hash_fq=Seq8_Hash_ATTACTCG_
#hash_fq=Seq13_Hash_ATCACGTT_
hash_fq=Seq18_Hash_ATCACGTT_
#hash_fq=Seq19_Hash_ATTACTAC_
#hash_fq=Seq20_Hash_ATCACGTT_,Seq20_Hash_ATTATACA_
#hash_fq=Seq21_Hash_ATTACTAC_
#hash_fq=Seq24_Hash_
#hash_fq=Seq25_Hash_
#hash_fq=Seq27_Hash_ATACACTG_
#hash_fq=Seq29_Hash_ATTATACA_
#hash_fq=Seq30_Hash_ATTACTAC_
#hash_fq=Seq31_Hash_ATCACGTT_

### Name of ECB masterlist to use
ecb_master=ECB_groups_master_list   # seq8, 18, 20
#ecb_master=ECB_groups_master_list_updated   # seq21, 24, 25

# multiple sample: separate by ',' with no space
# ecb_fq=Seq20_ECB_ATTACTCG_,Seq20_ECB_TCCGGAGA_
# hash_fq=Seq20_Hash_ATTATACA_,Seq20_Hash_ATCACGTT_


## To run (require python 3.3 or above) ##
pythonpath="/N/users/kasperk/anaconda3/bin/python"

$pythonpath ecb_parsing.py ${ecb_fqpath} ${ecb_fq} ${ecb_master}                                     ### Run ECB parsing before combined ECB and Hash analysis
$pythonpath hash_parsing.py ${hash_fqpath} ${hash_fq}   			        ### Run Hash parsing before combined ECB and Hash analysis. Run only this if no ECBs in sample

### RUN FOR SEQ8, 18, 20
$pythonpath combine_ecb_hash_bcrg_updated.py ${ecb_out} ${hash_out} ${ecb_fqpath}       ### Run for samples where the ECB masterlist should output updated version (i.e. RG combined based on dual insertions). SEQ8, 18, 20

### RUN FOR SEQ21, 24, 25
#$pythonpath combine_ecb_hash_bcrg_not_updated.py ${ecb_out} ${hash_out} ${ecb_fqpath}   ### Run for samples where the ECB masterlist should output non-updated version (i.e. ECB masterlist has been previously updated). SEQ21, 24, 25

