module rest

module load BWA/0.7.17-gimkl-2017a
module load SAMtools/1.9-GCC-7.4.0
module load BCFtools/1.9-GCC-7.4.0
module load BBMap/38.81-gimkl-2020a
module load R/4.0.1-gimkl-2020a
module load File-Rename/1.13-GCC-9.2.0


chmod -R +x /nesi/nobackup/massey02696/Azorella/Bioinformatics/HybPhaser/
chmod -R +x /scale_wlg_nobackup/filesets/nobackup/massey02696/Azorella/Bioinformatics/HybPhaser-main/
export PATH=$PATH:/nesi/nobackup/massey02696/Azorella/Bioinformatics/HybPhaser-main/:$PATH

/nesi/nobackup/massey02696/Azorella/Bioinformatics/HybPhaser-main/1_generate_consensus_sequences.sh \
-n /nesi/nobackup/massey02696/Azorella/Targetenrich/readsdata/HybPiper/namelist3.txt \
-p /nesi/nobackup/massey02696/Azorella/Targetenrich/readsdata/HybPiper \
-o /nesi/nobackup/massey02696/Azorella/Targetenrich/readsdata/HybPhaser -t 4 -i intronerated

R

library("ape")
library("seqinr")
library("stringr")


config_file="/nesi/nobackup/massey02696/Azorella/Targetenrich/readsdata/HybPhaser/config.txt"

##################################
### Part 1: Assessment of SNPs ###
##################################


# 1a) execute script to count SNPs in consensus files (this will take a few minutes)
source("/scale_wlg_nobackup/filesets/nobackup/massey02696/Azorella/Bioinformatics/HybPhaser-main/1a_count_snps.R")

# 1b) excute the script to generate graphs for the data assessment 
source("/scale_wlg_nobackup/filesets/nobackup/massey02696/Azorella/Bioinformatics/HybPhaser-main/1b_assess_dataset.R")

# 1c) generate sequence lists 
source("/scale_wlg_nobackup/filesets/nobackup/massey02696/Azorella/Bioinformatics/HybPhaser-main/1c_generate_sequence_lists.R")

# Sequence lists are available in subfolder and ready for alignment and phylogeny reconstruction
