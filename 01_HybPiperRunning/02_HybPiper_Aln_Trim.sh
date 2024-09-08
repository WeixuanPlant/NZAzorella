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

############################################################################################################

ls -d *nomerge > namelist.txt
hybpiper stats --targetfile_dna mega353.fasta supercontig namelist.txt
hybpiper paralog_retriever --targetfile_dna mega353.fasta namelist.txt
hybpiper recovery_heatmap seq_lengths.tsv

mkdir SupercontigsOutput
hybpiper retrieve_sequences --targetfile_dna mega353.fasta supercontig --sample_names namelist.txt --fasta_dir SupercontigsOutput

# Get a list of genes with sequences more than 15 indivdiuals 
grep ">" -c *_supercontig.fasta | sed 's/_supercontig.fasta:/ /' | sort -g -k 2 | awk '$2 >= 15  { print $1}' > aln132_n15.txt

# Remove all "NNNN" gap sites among sequences
for i in $(cat aln132_n15.txt); do  awk '{ print $1}' "${i}"_supercontig.fasta | \
 awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' | \
 sed '/^>/! s/N//g' > "${i}"_supercontig_fixed.fasta; done

# build a Tree rename file
#https://www.cyberciti.biz/faq/how-to-remove-carriage-return-in-linux-or-unix/
sed 's/\r$//' Treenamen132.csv > Treenamen132_fixed.csv


#######################################################################################

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=10G 
#SBATCH --time=02:00:00
#SBATCH --mail-user=EMAILHere
#SBATCH --mail-type=ALL
#SBATCH --open-mode=append
#SBATCH --output="aln.%J.out"
#SBATCH --job-name="aln"
#SBATCH --array=1-351%2

word=$(cat sample_n131_genelistn351.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

module load mafft/7.508

echo "${word}"

mafft --auto --thread 10 "${word}"_supercontig_n131.fasta > "${word}"_aln.fasta
/work/LAS/jfw-lab/weixuan/94_weixuan/n23_newanalysis/trimal/source/trimal -in "${word}"_aln.fasta -out "${word}"_trim.fasta -gt 0.7


grep -c '>' *_trim.fasta | sed 's/_trim.fasta:/ /' | sort -g -k 2 | awk '$2 >= 15  { print $1, $2}' > aln132_n15.txt

for word in $(cat ../aln132_n15.txt)
do
echo "${word}"
sed 's/,/\t/g' ../Treenamen132_fixed.csv | while read a b; do sed -i "s/$a/$b/" "${word}"_trim.fasta; done
done
	
