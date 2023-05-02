#!/bin/bash
## run the Rscript PANDA.R and schedule this job to SLURM with
## `sbatch PANDA.sh`

#SBATCH --job-name=PANDA
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jwhitlock@uab.edu #CHANGE THE EMAIL TO YOUR OWN
#
#Number of tasks needed for this job, the partition, and parameters
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=255000
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --partition=largemem #partition info here: https://docs.rc.uab.edu/cheaha/hardware/#partitions
#
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=0-49 #50 cell type inputs total across both conditions, 16 cortex and 34 kidney



########################################
### PUT YOUR COMMANDS BELOW THIS BOX ###
########################################

#load modules
module load Singularity/3.5.2-GCC-5.4.0-2.26

#variables
wd="$USER_DATA/230227_JW_Setbp1Manuscript"
src="$USER_DATA/230227_JW_Setbp1Manuscript/src/network_scripts/PANDA" #be sure that your subdirectories are structured the same

#code to execute docker and script for analysis
cd ${wd}

#array of cell type specific expression inputs for PANDA
SAMPLE_LIST="${wd}/results/array_inputs/Setbp1_PANDA_files.txt" #note: make sure path and file name are correct

singularity exec -B ${src} ${wd}/bin/PANDA_docker/setbp1_manuscript_panda_1.0.0_latest.sif Rscript --vanilla ${src}/PANDA.R ${SAMPLE_LIST[$SLURM_ARRAY_TASK_ID]} # here vanilla ensures only the script is run and environment is kept clean