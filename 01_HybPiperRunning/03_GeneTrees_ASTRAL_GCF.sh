#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G 
#SBATCH --time=0-05:00:00
#SBATCH --mail-type=ALL
#SBATCH --open-mode=append
#SBATCH --output="job.vcf_n319.%J.out"
#SBATCH --job-name="jointGeno_n319"
#SBATCH --array=1-351%2

module load iqtree2/2.2.2.7-py310-kcqit6n
file=$(ls -1 ../aln/*_trim.fasta | sed -n ${SLURM_ARRAY_TASK_ID}p)
name=$(basename $file _trim.fasta)
iqtree2 -s ../aln/${name}_trim.fasta --prefix $name -T AUTO -B 1000 -m MFP


# collapes the node bootstrapping value less than 10%
./nw_ed ./loci_gt07_n20_Azorella.treefile  'i & b<=10' o > loci_gt07_n20_Azorella.bs10.treefile

java -jar /scale_wlg_nobackup/filesets/nobackup/massey02696/Azorella/Bioinformatics/Astral/astral.5.7.7.jar -i loci_gt07_n20_Azorella.bs10.treefile -o loci_gt07_n20_Azorella.bs10.astral.treefile

# calculating the condordance factors
module load ETE/3.1.1-gimkl-2020a-Python-3.8.2

#make an outgroup.txt for rerooting trees
python3 reroot_trees.py loci_gt07_n20_Azorella.bs10.treefile outgroup.txt > loci_gt07_n20_Azorella.bs10.rerooted.treefile
python3 reroot_trees.py loci_gt07_n20_Azorella.bs10.astral.treefile outgroup.txt > loci_gt07_n20_Azorella.bs10.astral.rerooted.treefile
java -jar  -Xmx16g phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar -a 1 -v -d loci_gt07_n20_Azorella.bs10.rerooted.treefile -m loci_gt07_n20_Azorella.bs10.astral.rerooted.treefile  -o out_n20
 
grep 'Azly' loci_gt07_n4_Azorella.bs10.treefile  > loci_gt07_n4_Azorella.bs10.outgroup.treefile
python3 ../phypartspiecharts/reroot_trees.py loci_gt07_n4_Azorella.bs10.outgroup.treefile utgroup.list  > loci_gt07_n4_Azorella.bs10.outgroup.rerooted.treefile

#########################################################################################

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200G 
#SBATCH --time=3-05:00:00
#SBATCH --mail-type=ALL
#SBATCH --open-mode=append
#SBATCH --output="job.vcf_n319.%J.out"
#SBATCH --job-name="jointGeno_n319"

module load openjdk/17.0.5_8-5pdrwz2

java -jar -Xmx100g phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar -a 1 -v -d loci_gt07_n15_gene345_bs10_reroot.treefile -m loci_gt07_n15_gene345_bs10_reroot_astral_reroot.treefile -o Azorella131_n15_gene345_bs10













