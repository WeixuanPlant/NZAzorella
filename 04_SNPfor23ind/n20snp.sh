#!/bin/bash -e
#SBATCH --job-name=plastome_01  #job name (shows up in the queue)
#SBATCH --account=massey02696
#SBATCH --cpus-per-task=4   # number of CPUs per task (1 by default)
#SBATCH --time=24:00:00 #Walltime (HH:MM:SS)
#SBATCH --mem=100G  #Memory in GB 

module load Python/3.9.9-gimkl-2020a
module load BWA/0.7.17-gimkl-2017a
module load GATK/4.2.6.1-gimkl-2020a
module load PLINK/2.00a2.3
module load SAMtools/1.9-GCC-7.4.0
module load BCFtools/1.9-GCC-7.4.0

for i in $(cat namelist_n20_fixed.txt); do  echo "${i}"; bash variantcall.sh ref_wx5.fasta "${i}"; done

