#!/bin/bash -l

####################################
#     ARIS slurm script template   #
#                                  #
# Submit script: sbatch filename   #
#                                  #
####################################

#SBATCH --job-name=LASSO    # Job name
#SBATCH --output=/new-stg/home/banghua/TWAS_ASSOC/LASSO.%j.out # Stdout (%j expands to jobId)
#SBATCH --error=/new-stg/home/banghua/TWAS_ASSOC/LASSO.%j.err # Stderr (%j expands to jobId)
#SBATCH --ntasks=16     # Number of tasks(processes)
#SBATCH --nodes=1     # Number of nodes requested
#SBATCH --mem=128G   # memory per NODE
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=bax001@ucsd.edu # Email to which notifications will be sent

#write your own codes below:
conda activate CSE_284
cd /new-stg/home/banghua/TWAS_ASSOC
python /new-stg/home/banghua/TWAS_ASSOC/LASSO.py
