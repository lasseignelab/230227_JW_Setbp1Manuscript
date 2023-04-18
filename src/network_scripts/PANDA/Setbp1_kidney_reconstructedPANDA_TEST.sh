#!/bin/bash
#
#SBATCH --job-name=All_kid
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jwhitlock@uab.edu #CHANGE THE EMAIL TO YOUR OWN
#
#Number of tasks needed for this job, the partition, and parameters
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=255000
#SBATCH --nodes=1
#SBATCH --time=150:00:00
#SBATCH --partition=long
#
#SBATCH --output=Test23_kidney.out
#SBATCH --error=Test23_kidney.err


########################################
### PUT YOUR COMMANDS BELOW THIS BOX ###
########################################

#Set your environment here
module load R/4.2.0-foss-2021a

#script to call on
Rscript Setbp1_kidney_reconstructedPANDA_TEST.R
