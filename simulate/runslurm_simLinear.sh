#!/bin/bash -l
# call with:
# >> sbatch runslurm [submit as job]

# Job name
#SBATCH --job-name=30trialsPerCdt 
## Mail events (NONE, BEGIN, END, FAIL, ALL)
###############################################
########## example #SBATCH --mail-type=END,FAIL 
##############################################
##SBATCH --mail-type=FAIL
##SBATCH --mail-user=olg24@pitt.edu

# Run on a single CPU
#SBATCH --ntasks=1    
# Submit job to cpu queue                
#SBATCH -p cpu
#SBATCH -a 1-60%60

# Job memory request
#SBATCH --mem=10gb
# Time limit days-hrs:min:sec
#SBATCH --time 00-05:30:00

# Standard output and error log
#SBATCH -o ./log/matlabsim%A_%a.out

hostname
echo "job starting"
module load matlab-9.5
echo "RUNNING MATLAB"
matlab < Sim_1layer_cluster_sameParams.m
matlab -r "dumb; exit"
module unload matlab-9.5
echo "job finished"
