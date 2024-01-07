First, the sequencing reads were mapped to the Angiosperms353 bait set reference sequences (Johnson, et al. 2018) using BWA v. 0.7.17 (Li 2013), 
and to improve the exon recovery rate,

Second, mapped reads were sorted to each targeted gene and de novo assembled into large contigs with SPAdes v.3.11.1 (Bankevich, et al. 2012). 
Finally, the exons in assembled contigs were extracted and concatenated according to the reference sequences with Exonerate v.2.2.0 (Slater and Birney 2005). 
We also calculated the recovery rates of exons by comparing the recovered exon length with the length of references for each gene. 


#!/bin/bash -e
#SBATCH --job-name=angiosperms353_hybpiper  #job name (shows up in the queue)
#SBATCH --account=massey02696
#SBATCH --cpus-per-task=50   # number of CPUs per task (1 by default)
#SBATCH --time=100:00:00 #Walltime (HH:MM:SS)
#SBATCH --mem=50G  #Memory in GB 


module load Python/3.9.9-gimkl-2020a
module load HybPiper/2.0.1rc-Miniconda3
module load Parallel/20200522

for file in /nesi/nobackup/massey02696/WeixuanData/Azorella_Angiosperm353_NCBI_upload/02_pairedtrimed/*_1P.fq.gz
	do
	withpath="${file}"
	filename=${withpath##*/}
	base="${filename%*_1P.fq.gz}" 
	echo "${base}"
	
	hybpiper assemble --run_intronerate \
	--readfiles /nesi/nobackup/massey02696/WeixuanData/Azorella_Angiosperm353_NCBI_upload/02_pairedtrimed/"${base}"*P.fq.gz \
	--targetfile_dna mega353.fasta --bwa \
	--prefix "${base}"_nomerge  --no_padding_supercontigs \
	--timeout_assemble 600 --paralog_min_length_percentage 0.5
	done
	
