#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200G 
#SBATCH --time=1-05:00:00
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --mail-type=ALL
#SBATCH --open-mode=append
#SBATCH --output="job.vcf_n319.%J.out"
#SBATCH --job-name="jointGeno_n319"

module load openjdk/17.0.5_8-5pdrwz2

java -jar -Xmx100g phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar -a 1 -v -d loci_n219_bs10_reroot.treefile -m loci_n219_bs10_reroot_astral_reroot.treefile -o Azorella131_n23_gene219_bs10

