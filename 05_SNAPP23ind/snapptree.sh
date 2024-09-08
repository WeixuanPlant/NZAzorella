#!/bin/bash -e
#SBATCH --job-name=angiosperms353_snapp  #job name (shows up in the queue)
#SBATCH --account=massey02696
#SBATCH --cpus-per-task=35
#SBATCH --mem=10G  #Memory in GB 
#SBATCH --time=200:15:00      # Walltime (HH:MM:SS)

module load  BEAST/2.6.6
module load beagle-lib/4.0.0-GCC-11.3.0

beast -resume -instances 3 -threads 30 AZn23_SNAPP.length.xml
