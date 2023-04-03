#!/bin/bash

# Set job time to 3 days.
#SBATCH --time=3-00:00:00

# Set a name for the job, visible in `squeue`
#SBATCH --job-name=MY_JOB

# One node.
#SBATCH --nodes=1

# One task
#SBATCH --ntasks=1

# One CPU/core per task, n of threads
#SBATCH --cpus-per-task=8

# 4GB of RAM
#SBATCH --mem=64G

# Who to send mail to.
#SBATCH --mail-user=mjprzy@stanford.edu

# What type of mail to send
#SBATCH --mail-type=FAIL

# which account - this is essential
#SBATCH --account=ccurtis2

# load required modules which are need by your script
source activate snakes

# command to execute for n in $(seq 1 20); do sbatch /home/mjprzy/hdp_multichain.sh $n; done
outdir="/labs/ccurtis2/mjprzy/WGS/HDP"
n=$1
Rscript /home/mjprzy/hdp_single_chain.R $outdir $n