#!/bin/bash -e
#SBATCH --job-name=angiosperms353_aln  #job name (shows up in the queue)
#SBATCH --account=massey02696
#SBATCH --cpus-per-task=10   # number of CPUs per task (1 by default)
#SBATCH --time=10:00:00 #Walltime (HH:MM:SS)
#SBATCH --mem=10G  #Memory in GB 

module load IQ-TREE/2.2.0.5-gimpi-2022a
iqtree2 -S aln/ --prefix ./loci_gt07_n15_Azorella -T AUTO -B 1000 --undo


# collapes the node bootstrapping value less than 30%
./nw_ed ./loci_gt07_n20_Azorella.treefile  'i & b<=30' o > loci_gt07_n20_Azorella.bs30.treefile

java -jar /scale_wlg_nobackup/filesets/nobackup/massey02696/Azorella/Bioinformatics/Astral/astral.5.7.7.jar -i loci_gt07_n20_Azorella.bs30.treefile -o loci_gt07_n20_Azorella.bs30.astral.treefile

# calculating the condordance factors
module load ETE/3.1.1-gimkl-2020a-Python-3.8.2

#make an outgroup.txt for rerooting trees
python3 reroot_trees.py loci_gt07_n20_Azorella.bs30.treefile outgroup.txt > loci_gt07_n20_Azorella.bs30.rerooted.treefile
python3 reroot_trees.py loci_gt07_n20_Azorella.bs30.astral.treefile outgroup.txt > loci_gt07_n20_Azorella.bs30.astral.rerooted.treefile
java -jar  -Xmx16g phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar -a 1 -v -d loci_gt07_n20_Azorella.bs30.rerooted.treefile -m loci_gt07_n20_Azorella.bs30.astral.rerooted.treefile  -o out_n20
 
grep 'Azly' loci_gt07_n4_Azorella.bs30.treefile  > loci_gt07_n4_Azorella.bs30.outgroup.treefile
python3 ../phypartspiecharts/reroot_trees.py loci_gt07_n4_Azorella.bs30.outgroup.treefile utgroup.list  > loci_gt07_n4_Azorella.bs30.outgroup.rerooted.treefile














