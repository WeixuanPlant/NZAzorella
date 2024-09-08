#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G 
#SBATCH --time=0-02:00:00
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --mail-type=ALL
#SBATCH --open-mode=append
#SBATCH --output="job.vcf_n319.%J.out"
#SBATCH --job-name="jointGeno_n319"
#SBATCH --array=1-220%2

module load iqtree2/2.2.2.7-py310-kcqit6n
file=$(ls -1 ../aln/*_trim.fasta | sed -n ${SLURM_ARRAY_TASK_ID}p)
name=$(basename $file _trim.fasta)
iqtree2 -s ../aln/${name}_trim.fasta --prefix $name -T AUTO -B 1000 -m MFP
